#!/usr/bin/env python3
"""CMGDB benchmark harness.

Runs timed, correctness-validated benchmark scenarios covering the main
CMGDB usage modes (Python box maps, data-driven box maps, the C++ interval
map, and Conley index computation) on Leslie population models in 2D, 3D,
and 4D.

Every scenario carries a frozen reference of its mathematical output (Morse
sets, reachability edges, Conley indices). A run that does not reproduce its
reference exactly is reported as FAILED and its timings are not to be
trusted. References live in benchmarks/references.json and are regenerated
only with --update-refs, which is legitimate only when the mathematical
output is *intended* to change.

Usage:
  python benchmark.py                 # quick suite
  python benchmark.py --heavy        # heavy suite
  python benchmark.py --all          # both suites
  python benchmark.py --scenario leslie2d_python --repeat 5
  python benchmark.py --json out.json
  python benchmark.py --compare old.json new.json
  python benchmark.py --update-refs  # regenerate references for selected scenarios
  python benchmark.py --list

Each scenario runs in its own subprocess: this isolates peak-memory
measurement, silences the C++ core's stdout chatter, and keeps one crash
from killing the suite. Timings measure a fresh Model + compute call per
repetition. The Python map function is wrapped with a counter so results
split wall time into "inside the Python map" vs "C++ core" — the key
diagnostic for evaluating map-evaluation optimizations.
"""

import argparse
import gc
import json
import math
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_DIR = BENCH_DIR.parent
EXAMPLES_DIR = REPO_DIR / "examples"
DEFAULT_REFS_PATH = BENCH_DIR / "references.json"

# ---------------------------------------------------------------------------
# Leslie maps
#
# The d-dimensional Leslie population model
#   f(x)_0 = (sum_i theta_i x_i) * exp(-0.1 * sum_i x_i)
#   f(x)_k = 0.7 * x_{k-1}          for k = 1, ..., d-1
# with the classic CMGDB example parameters theta = (19.6, 23.68, ...).
# Domain bounds: f(x)_0 <= max_i(theta_i) * s * exp(-0.1 s) <= 23.68*10/e ~ 87.1
# and each further coordinate is contracted by 0.7, so the boxes below are
# forward-invariant overestimates of the attracting region.
# ---------------------------------------------------------------------------

LESLIE_THETA = [19.6, 23.68, 19.6, 23.68]
LESLIE_DOMAINS = {
    2: ([-0.001, -0.001], [90.0, 70.0]),
    3: ([-0.001, -0.001, -0.001], [90.0, 65.0, 45.0]),
    4: ([-0.001, -0.001, -0.001, -0.001], [90.0, 65.0, 45.0, 32.0]),
}


def leslie_point_map(dim):
    theta = LESLIE_THETA[:dim]

    def f(x):
        s = sum(x)
        y0 = sum(t * xi for t, xi in zip(theta, x)) * math.exp(-0.1 * s)
        return [y0] + [0.7 * xi for xi in x[:dim - 1]]

    return f


def leslie_point_map_vectorized(dim):
    """Vectorized Leslie map for generating synthetic datasets: X (n, dim) -> (n, dim)."""
    import numpy as np
    theta = np.array(LESLIE_THETA[:dim])

    def f(X):
        s = X.sum(axis=1)
        y0 = (X @ theta) * np.exp(-0.1 * s)
        return np.column_stack([y0] + [0.7 * X[:, i] for i in range(dim - 1)])

    return f


# ---------------------------------------------------------------------------
# Map-call instrumentation
# ---------------------------------------------------------------------------

class CountingMap:
    """Wraps a rect -> rect box map, recording call count and cumulative time."""

    def __init__(self, F):
        self._F = F
        self.calls = 0
        self.seconds = 0.0

    def __call__(self, rect):
        t0 = time.perf_counter()
        result = self._F(rect)
        self.seconds += time.perf_counter() - t0
        self.calls += 1
        return result


class CountingBatchMap:
    """Wraps a batched box map, recording rectangles processed and time.

    calls counts rectangles (not invocations) so the "map calls" column is
    comparable between scalar and batched scenarios."""

    def __init__(self, F_batch):
        self._F = F_batch
        self.calls = 0
        self.seconds = 0.0

    def __call__(self, rects):
        t0 = time.perf_counter()
        result = self._F(rects)
        self.seconds += time.perf_counter() - t0
        self.calls += len(rects)
        return result


# ---------------------------------------------------------------------------
# Reference extraction
#
# The reference pins the mathematical output of a run in a form invariant
# under relabeling of Morse graph vertices: vertices are put in canonical
# order by the minimum box index of their Morse sets (box indices are leaf
# ranks of the grid tree, which are deterministic for a given grid), and all
# data is expressed in canonical labels. edges_unreduced is the full
# reachability relation -- the strongest invariant available.
# ---------------------------------------------------------------------------

def extract_reference(morse_graph, map_graph=None):
    import CMGDB  # noqa: F401  (morse_graph is a CMGDB object)
    num_vertices = morse_graph.num_vertices()
    morse_sets = [sorted(morse_graph.morse_set(v)) for v in range(num_vertices)]
    canonical = sorted(range(num_vertices), key=lambda v: morse_sets[v][0])
    relabel = {v: i for i, v in enumerate(canonical)}

    reference = {
        "num_vertices": num_vertices,
        "morse_set_sizes": [len(morse_sets[v]) for v in canonical],
        "morse_set_min_box": [morse_sets[v][0] for v in canonical],
        "edges_unreduced": sorted(
            [relabel[u], relabel[w]] for u, w in morse_graph.edges_unreduced()
        ),
    }
    annotations = [list(morse_graph.annotations(v)) for v in canonical]
    if any(annotations):
        reference["conley_indices"] = annotations
    if map_graph is not None:
        reference["phase_space_size"] = map_graph.num_vertices()
    return reference


# ---------------------------------------------------------------------------
# Scenario definitions
#
# Each scenario is a function returning (result_reference, metrics) for one
# full computation. Subdivision depths are calibrated so the quick suite
# runs in a couple of minutes total; heavy scenarios are meaningfully larger.
# ---------------------------------------------------------------------------

def _run_python_leslie(dim, subdiv_min, subdiv_max, conley=False, batch=False):
    import CMGDB
    f = leslie_point_map(dim)
    F = CountingMap(lambda rect: CMGDB.BoxMap(f, rect))
    lower_bounds, upper_bounds = LESLIE_DOMAINS[dim]
    model = CMGDB.Model(subdiv_min, subdiv_max, lower_bounds, upper_bounds, F)
    F_batch = None
    if batch:
        f_vec = leslie_point_map_vectorized(dim)
        F_batch = CountingBatchMap(lambda rects: CMGDB.BoxMapBatch(f_vec, rects))
        model.set_batch_map(F_batch)
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    reference = extract_reference(morse_graph, map_graph)
    map_calls = F.calls + (F_batch.calls if F_batch else 0)
    map_seconds = F.seconds + (F_batch.seconds if F_batch else 0.0)
    return reference, {"map_calls": map_calls, "map_seconds": map_seconds}


def _run_data_leslie(subdiv_min, subdiv_max, num_points=None, seed=42,
                     linear=False, batch=False, conley=False):
    """BoxMapData scenario. With num_points=None uses examples/Leslie_500.dat;
    otherwise generates a deterministic synthetic dataset of that size.
    linear selects the linear-scan reference implementation; batch attaches
    the indexed implementation's batch evaluator via set_batch_map."""
    import numpy as np
    import CMGDB
    lower_bounds, upper_bounds = LESLIE_DOMAINS[2]
    if num_points is None:
        data = np.loadtxt(EXAMPLES_DIR / "Leslie_500.dat")
        X, Y = data[:, :2], data[:, 2:]
    else:
        rng = np.random.default_rng(seed)
        X = rng.uniform(lower_bounds, upper_bounds, size=(num_points, 2))
        Y = leslie_point_map_vectorized(2)(X)
    data_map = (CMGDB.BoxMapDataLinear if linear else CMGDB.BoxMapData)(X, Y)
    F = CountingMap(data_map)
    model = CMGDB.Model(subdiv_min, subdiv_max, lower_bounds, upper_bounds, F)
    F_batch = None
    if batch:
        F_batch = CountingBatchMap(data_map.batch)
        model.set_batch_map(F_batch)
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    reference = extract_reference(morse_graph, map_graph)
    map_calls = F.calls + (F_batch.calls if F_batch else 0)
    map_seconds = F.seconds + (F_batch.seconds if F_batch else 0.0)
    return reference, {"map_calls": map_calls, "map_seconds": map_seconds}


def _make_expensive_leslie2d(width=1024, seed=7, amplitude=0.5):
    """Leslie map plus a small random-MLP perturbation: keeps Leslie-like
    dynamics (so the adaptive tree explores a substantial recurrent region)
    while making each evaluation cost two (m, width) matrix products --
    a stand-in for neural-network / Gaussian-process surrogate maps."""
    import numpy as np
    rng = np.random.default_rng(seed)
    W1 = rng.normal(0.0, 1.0, (2, width))
    b1 = rng.normal(0.0, 1.0, width)
    W2 = rng.normal(0.0, 1.0 / np.sqrt(width), (width, 2))
    b2 = rng.normal(0.0, 1.0, 2)
    scale_in = np.array([45.0, 35.0])
    base = leslie_point_map_vectorized(2)

    def f(P):
        hidden = np.tanh(P / scale_in @ W1 + b1)
        return base(P) + amplitude * np.tanh(hidden @ W2 + b2)

    return f


def _run_expensive_leslie2d(subdiv_min, subdiv_max, precomputed=False, width=1024):
    import numpy as np
    import CMGDB
    lower_bounds, upper_bounds = LESLIE_DOMAINS[2]
    f = _make_expensive_leslie2d(width)

    def f_scalar(x):
        return list(f(np.array([x]))[0])

    F = CountingMap(lambda rect: CMGDB.BoxMap(f_scalar, rect))
    model = CMGDB.Model(subdiv_min, subdiv_max, lower_bounds, upper_bounds, F)
    if precomputed:
        # The precompute cost is inside the timed region: this measures the
        # honest end-to-end trade (build table once, then lookups only).
        # batch_points bounds the hidden-layer temporaries of the wide MLP.
        F_pre = CMGDB.PrecomputedBoxMap(f, lower_bounds, upper_bounds, subdiv_max,
                                        batch_points=16384)
        F_batch = CountingBatchMap(F_pre.batch)
    else:
        F_batch = CountingBatchMap(lambda rects: CMGDB.BoxMapBatch(f, rects))
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model, batch_chunk_size=16384)
    reference = extract_reference(morse_graph, map_graph)
    map_calls = F.calls + F_batch.calls
    map_seconds = F.seconds + F_batch.seconds
    return reference, {"map_calls": map_calls, "map_seconds": map_seconds}


def _run_product(dim, subdiv_min, subdiv_max, subdiv_init, conley=False):
    """The decoupled product map f(x)_i = x_i/(2-x_i) on [0, 1.2]^dim
    (the tests/test_basic.py map). Each coordinate has attracting fixed points
    at 0 and 1, so the Morse graph is a product partial order with 5^dim
    Morse sets under uniform decomposition -- a richer edge structure than
    the Leslie attractor scenarios (and in 3D more than 64 Morse sets, which
    exercises the multi-word reachability grouping). With
    subdiv_min == subdiv_max == subdiv_init this exercises the uniform
    (fixed-depth) decomposition path."""
    import CMGDB

    def f(x):
        return [x[d] / (2.0 - x[d]) for d in range(dim)]

    F = CountingMap(lambda rect: CMGDB.BoxMap(f, rect))
    model = CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, 10000,
                        [0.0] * dim, [1.2] * dim, F)
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    reference = extract_reference(morse_graph, map_graph)
    return reference, {"map_calls": F.calls, "map_seconds": F.seconds}


def _run_henon3d(subdiv_min, subdiv_max, conley=False):
    """Delay-embedded Henon map in 3D: a chaotic attractor whose recurrent
    set stays large under subdivision -- the opposite stress profile from the
    strongly contracting Leslie scenarios (big SCCs, heavy cover load, and a
    large index pair for the Conley variant's 3D homology)."""
    import CMGDB

    def f(x):
        return [1.0 - 1.4 * x[0] * x[0] + 0.3 * x[1], x[0], x[1]]

    # padding=True: corner sampling alone cannot bound the quadratic fold
    # (its maximum is interior to the box), so pad to an outer enclosure
    F = CountingMap(lambda rect: CMGDB.BoxMap(f, rect, padding=True))
    model = CMGDB.Model(subdiv_min, subdiv_max,
                        [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], F)
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    reference = extract_reference(morse_graph, map_graph)
    return reference, {"map_calls": F.calls, "map_seconds": F.seconds}


def _run_interval_leslie(subdiv_min, subdiv_max):
    """The built-in C++ interval-arithmetic Leslie map: no Python callback."""
    import CMGDB
    lower_bounds, upper_bounds = LESLIE_DOMAINS[2]
    morse_graph = CMGDB.MorseGraphIntvalMap(
        subdiv_min, subdiv_max, lower_bounds, upper_bounds,
        [19.6, 23.68], "morse_intval.csv")
    reference = extract_reference(morse_graph)
    return reference, {"map_calls": 0, "map_seconds": 0.0}


SCENARIOS = {
    # ---- quick suite ----
    "leslie2d_python": {
        "suite": "quick",
        "run": lambda: _run_python_leslie(2, 20, 24),
        "describe": "2D Leslie, Python map via BoxMap corners, subdiv 20/24",
    },
    "leslie3d_python": {
        "suite": "quick",
        "run": lambda: _run_python_leslie(3, 21, 24),
        "describe": "3D Leslie, Python map via BoxMap corners, subdiv 21/24",
    },
    "leslie2d_data": {
        "suite": "quick",
        "run": lambda: _run_data_leslie(15, 18),
        "describe": "2D Leslie, BoxMapData with Leslie_500.dat, subdiv 15/18",
    },
    "leslie2d_data_linear": {
        "suite": "quick",
        "run": lambda: _run_data_leslie(15, 18, linear=True),
        "describe": "2D Leslie, linear-scan BoxMapDataLinear with Leslie_500.dat, subdiv 15/18",
    },
    "leslie2d_interval": {
        "suite": "quick",
        "run": lambda: _run_interval_leslie(20, 26),
        "describe": "2D Leslie, built-in C++ interval map, subdiv 20/26",
    },
    "leslie2d_conley": {
        "suite": "quick",
        "run": lambda: _run_python_leslie(2, 20, 22, conley=True),
        "describe": "2D Leslie, Conley-Morse graph, subdiv 20/22",
    },
    "product2d_adaptive": {
        "suite": "quick",
        "run": lambda: _run_product(2, 14, 18, 8),
        "describe": "2D product map x/(2-x), adaptive subdiv 14/18 with init 8, product-order Morse graph",
    },
    "product2d_uniform": {
        "suite": "quick",
        "run": lambda: _run_product(2, 12, 12, 12),
        "describe": "2D product map, uniform (fixed-depth) decomposition at subdiv 12, 25 Morse sets",
    },
    "product2d_conley": {
        "suite": "quick",
        "run": lambda: _run_product(2, 14, 18, 8, conley=True),
        "describe": "2D product map, Conley-Morse graph, subdiv 14/18 init 8",
    },
    "product3d_uniform": {
        "suite": "quick",
        "run": lambda: _run_product(3, 15, 15, 15),
        "describe": "3D product map, uniform decomposition at subdiv 15, 125 Morse sets (>64: multi-word reachability)",
    },
    "product4d_adaptive": {
        "suite": "quick",
        "run": lambda: _run_product(4, 16, 18, 14),
        "describe": "4D product map, adaptive subdiv 16/18 with init 14, 225 Morse sets (>64 under adaptive decomposition)",
    },
    "henon3d_conley": {
        "suite": "quick",
        "run": lambda: _run_henon3d(12, 15, conley=True),
        "describe": "3D Henon attractor, Conley-Morse graph (3D homology), subdiv 12/15",
    },
    "leslie2d_data_conley": {
        "suite": "quick",
        "run": lambda: _run_data_leslie(15, 18, conley=True),
        "describe": "2D Leslie, BoxMapData with Leslie_500.dat, Conley-Morse graph, subdiv 15/18",
    },
    "leslie2d_python_batch": {
        "suite": "quick",
        "run": lambda: _run_python_leslie(2, 20, 24, batch=True),
        "describe": "2D Leslie, batched NumPy map via set_batch_map, subdiv 20/24",
    },
    "leslie3d_python_batch": {
        "suite": "quick",
        "run": lambda: _run_python_leslie(3, 21, 24, batch=True),
        "describe": "3D Leslie, batched NumPy map via set_batch_map, subdiv 21/24",
    },
    # ---- heavy suite ----
    "leslie2d_python_deep": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(2, 24, 28),
        "describe": "2D Leslie, Python map, subdiv 24/28",
    },
    "leslie3d_python_deep": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(3, 24, 27),
        "describe": "3D Leslie, Python map, subdiv 24/27",
    },
    "leslie4d_python": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(4, 22, 26),
        "describe": "4D Leslie, Python map, subdiv 22/26",
    },
    "leslie2d_data_large": {
        "suite": "heavy",
        "run": lambda: _run_data_leslie(14, 16, num_points=50_000),
        "describe": "2D Leslie, BoxMapData with 50k synthetic points, subdiv 14/16",
    },
    "leslie2d_python_deep_batch": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(2, 24, 28, batch=True),
        "describe": "2D Leslie, batched NumPy map, subdiv 24/28",
    },
    "leslie2d_data_large_batch": {
        "suite": "heavy",
        "run": lambda: _run_data_leslie(14, 16, num_points=50_000, batch=True),
        "describe": "2D Leslie, indexed BoxMapData batch path, 50k points, subdiv 14/16",
    },
    "henon3d": {
        "suite": "heavy",
        "run": lambda: _run_henon3d(21, 24),
        "describe": "3D Henon attractor (large chaotic recurrent set), subdiv 21/24",
    },
    "product3d_uniform_deep": {
        "suite": "heavy",
        "run": lambda: _run_product(3, 18, 18, 18),
        "describe": "3D product map, uniform decomposition at subdiv 18, 125 Morse sets over 262k boxes",
    },
    "leslie3d_conley": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(3, 18, 21, conley=True),
        "describe": "3D Leslie, Conley-Morse graph (3D homology), subdiv 18/21",
    },
    "leslie4d_conley": {
        "suite": "heavy",
        "run": lambda: _run_python_leslie(4, 16, 20, conley=True),
        "describe": "4D Leslie, Conley-Morse graph (4D homology), subdiv 16/20",
    },
    "leslie2d_expensive_live": {
        "suite": "heavy",
        "run": lambda: _run_expensive_leslie2d(15, 16),
        "describe": "2D Leslie + MLP perturbation (expensive map), live batched evaluation, subdiv 15/16",
    },
    "leslie2d_expensive_pre": {
        "suite": "heavy",
        "run": lambda: _run_expensive_leslie2d(15, 16, precomputed=True),
        "describe": "2D Leslie + MLP perturbation (expensive map), PrecomputedBoxMap, subdiv 15/16",
    },
}

DEFAULT_REPEAT = {"quick": 3, "heavy": 1}


# ---------------------------------------------------------------------------
# Child process: run one scenario, write JSON
# ---------------------------------------------------------------------------

def run_child(name, repeat, out_path):
    scenario = SCENARIOS[name]
    wall_times = []
    map_calls = []
    map_seconds = []
    reference = None
    for _ in range(repeat):
        gc.collect()
        t0 = time.perf_counter()
        ref, metrics = scenario["run"]()
        wall_times.append(time.perf_counter() - t0)
        map_calls.append(metrics["map_calls"])
        map_seconds.append(metrics["map_seconds"])
        if reference is None:
            reference = ref
        elif ref != reference:
            print(f"ERROR: scenario {name} is nondeterministic across repetitions",
                  file=sys.stderr)
            return 2
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is bytes on macOS, kilobytes on Linux
    maxrss_mb = maxrss / 1e6 if sys.platform == "darwin" else maxrss / 1e3
    payload = {
        "scenario": name,
        "wall_times": wall_times,
        "map_calls": map_calls,
        "map_seconds": map_seconds,
        "maxrss_mb": maxrss_mb,
        "reference": reference,
    }
    Path(out_path).write_text(json.dumps(payload))
    return 0


# ---------------------------------------------------------------------------
# Parent process: orchestrate scenarios, validate, report
# ---------------------------------------------------------------------------

def collect_metadata():
    meta = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
    }
    try:
        import numpy
        meta["numpy"] = numpy.__version__
    except ImportError:
        pass
    try:
        import importlib.metadata
        meta["cmgdb_version"] = importlib.metadata.version("CMGDB")
    except Exception:
        pass
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_DIR,
            capture_output=True, text=True, timeout=10)
        if sha.returncode == 0:
            meta["git_sha"] = sha.stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO_DIR,
            capture_output=True, text=True, timeout=10)
        if dirty.returncode == 0:
            meta["git_dirty"] = bool(dirty.stdout.strip())
    except Exception:
        pass
    if sys.platform == "darwin":
        try:
            cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                 capture_output=True, text=True, timeout=10)
            if cpu.returncode == 0:
                meta["cpu"] = cpu.stdout.strip()
        except Exception:
            pass
    return meta


def run_scenario_subprocess(name, repeat, timeout):
    with tempfile.TemporaryDirectory(prefix="cmgdb_bench_") as tmpdir:
        out_path = Path(tmpdir) / "result.json"
        cmd = [sys.executable, str(Path(__file__).resolve()),
               "--child", name, "--repeat", str(repeat), "--out", str(out_path)]
        try:
            proc = subprocess.run(
                cmd, cwd=tmpdir, stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return None, f"timed out after {timeout}s"
        if proc.returncode != 0 or not out_path.exists():
            stderr_tail = "\n".join(proc.stderr.strip().splitlines()[-8:])
            return None, f"exit code {proc.returncode}\n{stderr_tail}"
        return json.loads(out_path.read_text()), None


def diff_references(expected, actual):
    """Human-readable summary of how two references differ."""
    lines = []
    keys = sorted(set(expected) | set(actual))
    for key in keys:
        exp, act = expected.get(key), actual.get(key)
        if exp == act:
            continue
        if isinstance(exp, list) and isinstance(act, list) and len(exp) != len(act):
            lines.append(f"  {key}: length {len(exp)} -> {len(act)}")
        else:
            exp_s, act_s = str(exp), str(act)
            if len(exp_s) > 60:
                exp_s = exp_s[:60] + "..."
            if len(act_s) > 60:
                act_s = act_s[:60] + "..."
            lines.append(f"  {key}: {exp_s} -> {act_s}")
    return "\n".join(lines) or "  (no field differences found)"


def median(values):
    values = sorted(values)
    n = len(values)
    mid = n // 2
    return values[mid] if n % 2 else (values[mid - 1] + values[mid]) / 2


def fmt_time(seconds):
    if seconds < 10:
        return f"{seconds:.2f}s"
    if seconds < 120:
        return f"{seconds:.1f}s"
    return f"{seconds / 60:.1f}m"


def print_results_table(results):
    header = (f"{'scenario':<22} {'status':<8} {'wall (median)':<14} "
              f"{'map%':>5} {'map calls':>10} {'boxes':>10} {'sets':>5} {'MB':>7}")
    print(header)
    print("-" * len(header))
    for name, entry in results.items():
        if entry.get("error"):
            print(f"{name:<22} {'ERROR':<8} {entry['error'].splitlines()[0]}")
            continue
        status = entry["status"]
        wall = median(entry["wall_times"])
        map_sec = median(entry["map_seconds"])
        map_pct = 100.0 * map_sec / wall if wall > 0 else 0.0
        ref = entry["reference"]
        boxes = ref.get("phase_space_size", "-")
        print(f"{name:<22} {status:<8} {fmt_time(wall):<14} "
              f"{map_pct:>4.0f}% {max(entry['map_calls']):>10} "
              f"{boxes:>10} {ref['num_vertices']:>5} {entry['maxrss_mb']:>7.0f}")


def compare_files(path_a, path_b):
    a = json.loads(Path(path_a).read_text())
    b = json.loads(Path(path_b).read_text())
    header = (f"{'scenario':<22} {'old':<12} {'new':<12} {'speedup':>8}  "
              f"{'old map%':>8} {'new map%':>8}")
    print(f"old: {path_a} ({a['meta'].get('git_sha', '?')})")
    print(f"new: {path_b} ({b['meta'].get('git_sha', '?')})")
    print(header)
    print("-" * len(header))
    for name in a["results"]:
        if name not in b["results"]:
            continue
        ra, rb = a["results"][name], b["results"][name]
        if ra.get("error") or rb.get("error"):
            continue
        wa, wb = median(ra["wall_times"]), median(rb["wall_times"])
        pa = 100.0 * median(ra["map_seconds"]) / wa if wa else 0.0
        pb = 100.0 * median(rb["map_seconds"]) / wb if wb else 0.0
        flag = "" if ra["status"] == rb["status"] == "OK" else "  [CHECK STATUS]"
        print(f"{name:<22} {fmt_time(wa):<12} {fmt_time(wb):<12} "
              f"{wa / wb:>7.2f}x  {pa:>7.0f}% {pb:>7.0f}%{flag}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--heavy", action="store_true", help="run heavy suite")
    parser.add_argument("--all", action="store_true", help="run quick + heavy suites")
    parser.add_argument("--scenario", action="append", help="run specific scenario(s)")
    parser.add_argument("--repeat", type=int, help="repetitions per scenario")
    parser.add_argument("--json", metavar="PATH", help="write results to JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("OLD", "NEW"),
                        help="compare two results files and exit")
    parser.add_argument("--update-refs", action="store_true",
                        help="write references for the selected scenarios")
    parser.add_argument("--refs", default=str(DEFAULT_REFS_PATH),
                        help="path to references file")
    parser.add_argument("--no-validate", action="store_true",
                        help="skip reference validation (calibration only)")
    parser.add_argument("--timeout", type=int, default=3600,
                        help="per-scenario timeout in seconds")
    parser.add_argument("--list", action="store_true", help="list scenarios and exit")
    # child-mode internals
    parser.add_argument("--child", help=argparse.SUPPRESS)
    parser.add_argument("--out", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.child:
        sys.exit(run_child(args.child, args.repeat or 1, args.out))

    if args.compare:
        compare_files(*args.compare)
        return

    if args.list:
        for name, sc in SCENARIOS.items():
            print(f"{name:<22} [{sc['suite']}]  {sc['describe']}")
        return

    if args.scenario:
        unknown = [s for s in args.scenario if s not in SCENARIOS]
        if unknown:
            sys.exit(f"unknown scenario(s): {', '.join(unknown)}")
        selected = args.scenario
    elif args.all:
        selected = list(SCENARIOS)
    elif args.heavy:
        selected = [s for s, sc in SCENARIOS.items() if sc["suite"] == "heavy"]
    else:
        selected = [s for s, sc in SCENARIOS.items() if sc["suite"] == "quick"]

    refs_path = Path(args.refs)
    references = json.loads(refs_path.read_text()) if refs_path.exists() else {}

    results = {}
    failed = False
    for name in selected:
        suite = SCENARIOS[name]["suite"]
        repeat = args.repeat or DEFAULT_REPEAT[suite]
        print(f"running {name} (repeat={repeat}) ...", flush=True)
        payload, error = run_scenario_subprocess(name, repeat, args.timeout)
        if error:
            results[name] = {"error": error}
            failed = True
            print(f"  ERROR: {error}")
            continue
        if args.update_refs:
            references[name] = payload["reference"]
            payload["status"] = "NEW-REF"
        elif args.no_validate:
            payload["status"] = "NOVAL"
        elif name not in references:
            payload["status"] = "NO-REF"
            failed = True
            print(f"  no reference for {name}; run with --update-refs to create one")
        elif payload["reference"] != references[name]:
            payload["status"] = "FAILED"
            failed = True
            print(f"  VALIDATION FAILED for {name}:")
            print(diff_references(references[name], payload["reference"]))
        else:
            payload["status"] = "OK"
        results[name] = payload

    if args.update_refs:
        refs_path.write_text(json.dumps(references, indent=1))
        print(f"references written to {refs_path}")

    print()
    print_results_table(results)

    if args.json:
        output = {"meta": collect_metadata(), "results": results}
        Path(args.json).write_text(json.dumps(output, indent=1))
        print(f"\nresults written to {args.json}")

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
