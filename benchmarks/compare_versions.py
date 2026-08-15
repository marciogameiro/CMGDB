#!/usr/bin/env python3
"""Compare CMGDB performance between two versions of the code.

Builds the current working tree and an older git revision into two isolated
virtual environments, runs an identical scenario suite under both, and
writes a Markdown speedup table. Mechanisms that do not exist in the old
revision (batched maps, data-map batch, precomputed box maps) are
feature-detected and reported as n/a there.

If PyTorch is importable in the new version's environment, a third column
runs every scenario on the new version with the map evaluated by torch
(float64, so all columns compute identical Morse sets; float32 would be
faster still but can perturb box covers). Scenarios whose map cannot be
torch-evaluated (the C++ interval map; data-defined maps, which evaluate
no function) are n/a in that column, and scalar scenarios run with a
torch-backed batched map attached.

Usage:
  python benchmarks/compare_versions.py                  # vs 857ec8b (pre-optimization)
  python benchmarks/compare_versions.py --old-rev v1.3.2 # vs any git revision
  python benchmarks/compare_versions.py --scenarios leslie2d_python henon3d
  python benchmarks/compare_versions.py --no-torch       # skip the torch column
  python benchmarks/compare_versions.py --out my_table.md

The two builds take a few minutes; results land in version_compare.md (or
--out) next to this script. Scenarios run one repetition each, all flavors
back-to-back per scenario for thermal fairness. The scenario suite mirrors
benchmarks/benchmark.py but validates only that all flavors compute the
same number of Morse sets (the frozen references in references.json remain
the authoritative correctness gate for the current code).
"""

import argparse
import json
import math
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EXAMPLES = REPO / "examples"
DEFAULT_OLD_REV = "857ec8b"
TIMEOUT = 1800

USE_TORCH = False  # set in child mode by --use-torch

# ---------------------------------------------------------------------------
# Maps (defined here so both versions run identical code). Each map "kind"
# has a scalar form, a NumPy-vectorized form, and a torch float64 form.
# ---------------------------------------------------------------------------

LESLIE_THETA = [19.6, 23.68, 19.6, 23.68]
LESLIE_DOMAINS = {
    2: ([-0.001, -0.001], [90.0, 70.0]),
    3: ([-0.001, -0.001, -0.001], [90.0, 65.0, 45.0]),
    4: ([-0.001, -0.001, -0.001, -0.001], [90.0, 65.0, 45.0, 32.0]),
}


def scalar_map(kind, dim):
    if kind == "leslie":
        theta = LESLIE_THETA[:dim]

        def f(x):
            s = sum(x)
            y0 = sum(t * xi for t, xi in zip(theta, x)) * math.exp(-0.1 * s)
            return [y0] + [0.7 * xi for xi in x[:dim - 1]]
    elif kind == "product":
        def f(x):
            return [x[d] / (2.0 - x[d]) for d in range(dim)]
    elif kind == "henon":
        def f(x):
            return [1.0 - 1.4 * x[0] * x[0] + 0.3 * x[1], x[0], x[1]]
    else:
        raise ValueError(kind)
    return f


def numpy_map(kind, dim):
    import numpy as np
    if kind == "leslie":
        theta = np.array(LESLIE_THETA[:dim])

        def f(X):
            s = X.sum(axis=1)
            y0 = (X @ theta) * np.exp(-0.1 * s)
            return np.column_stack([y0] + [0.7 * X[:, i] for i in range(dim - 1)])
    elif kind == "product":
        def f(X):
            return X / (2.0 - X)
    elif kind == "henon":
        def f(X):
            return np.column_stack(
                [1.0 - 1.4 * X[:, 0] * X[:, 0] + 0.3 * X[:, 1], X[:, 0], X[:, 1]])
    else:
        raise ValueError(kind)
    return f


def torch_map(kind, dim):
    import torch

    if kind == "leslie":
        theta = torch.tensor(LESLIE_THETA[:dim], dtype=torch.float64)

        def f(X):
            with torch.no_grad():
                T = torch.as_tensor(X, dtype=torch.float64)
                s = T.sum(dim=1)
                y0 = (T @ theta) * torch.exp(-0.1 * s)
                cols = [y0] + [0.7 * T[:, i] for i in range(dim - 1)]
                return torch.stack(cols, dim=1).numpy()
    elif kind == "product":
        def f(X):
            with torch.no_grad():
                T = torch.as_tensor(X, dtype=torch.float64)
                return (T / (2.0 - T)).numpy()
    elif kind == "henon":
        def f(X):
            with torch.no_grad():
                T = torch.as_tensor(X, dtype=torch.float64)
                return torch.stack(
                    [1.0 - 1.4 * T[:, 0] * T[:, 0] + 0.3 * T[:, 1],
                     T[:, 0], T[:, 1]], dim=1).numpy()
    else:
        raise ValueError(kind)
    return f


def vec_map(kind, dim):
    return torch_map(kind, dim) if USE_TORCH else numpy_map(kind, dim)


def expensive_leslie2d_vec(width=1024, seed=7, amplitude=0.5):
    """Leslie + random-MLP perturbation: NumPy and torch float64 variants
    share the same weights, so all flavors compute the same map."""
    import numpy as np
    rng = np.random.default_rng(seed)
    W1 = rng.normal(0.0, 1.0, (2, width))
    b1 = rng.normal(0.0, 1.0, width)
    W2 = rng.normal(0.0, 1.0 / np.sqrt(width), (width, 2))
    b2 = rng.normal(0.0, 1.0, 2)
    scale_in = np.array([45.0, 35.0])
    base_np = numpy_map("leslie", 2)

    if USE_TORCH:
        import torch
        tW1 = torch.tensor(W1); tb1 = torch.tensor(b1)
        tW2 = torch.tensor(W2); tb2 = torch.tensor(b2)
        tscale = torch.tensor(scale_in)
        base_t = torch_map("leslie", 2)

        def f(P):
            with torch.no_grad():
                T = torch.as_tensor(P, dtype=torch.float64)
                hidden = torch.tanh(T / tscale @ tW1 + tb1)
                pert = amplitude * torch.tanh(hidden @ tW2 + tb2)
                return base_t(P) + pert.numpy()
        return f

    def f(P):
        hidden = np.tanh(P / scale_in @ W1 + b1)
        return base_np(P) + amplitude * np.tanh(hidden @ W2 + b2)

    return f


def box_map_batch(f_vec, rects, padding=False):
    import numpy as np
    rects = np.asarray(rects, dtype=float)
    N, two_dim = rects.shape
    dim = two_dim // 2
    lower, upper = rects[:, :dim], rects[:, dim:]
    num_corners = 2 ** dim
    corners = np.empty((num_corners, N, dim))
    for k in range(num_corners):
        mask = np.array([(k >> d) & 1 for d in range(dim)], dtype=bool)
        corners[k] = np.where(mask, upper, lower)
    Y = np.asarray(f_vec(corners.reshape(num_corners * N, dim)))
    Y = Y.reshape(num_corners, N, dim)
    Y_lower, Y_upper = Y.min(axis=0), Y.max(axis=0)
    if padding:
        pad = upper - lower
        Y_lower, Y_upper = Y_lower - pad, Y_upper + pad
    return np.hstack([Y_lower, Y_upper])


# ---------------------------------------------------------------------------
# Scenario implementations (feature-detecting)
# ---------------------------------------------------------------------------

class Skip(Exception):
    pass


def run_scalar(CMGDB, kind, dim, lower, upper, subdiv, conley=False,
               padding=False, batch=False):
    """Generic Model + ComputeMorseGraph scenario. In the torch flavor a
    torch-backed batched map is attached even for scalar scenarios (per-box
    torch calls would measure only call overhead)."""
    f = scalar_map(kind, dim)

    def F(rect):
        return CMGDB.BoxMap(f, rect, padding=padding)

    if len(subdiv) == 2:
        model = CMGDB.Model(subdiv[0], subdiv[1], lower, upper, F)
    else:
        model = CMGDB.Model(subdiv[0], subdiv[1], subdiv[2], 10000,
                            lower, upper, F)
    if batch or USE_TORCH:
        if not hasattr(model, "set_batch_map"):
            raise Skip("no set_batch_map")
        fv = vec_map(kind, dim)
        model.set_batch_map(lambda rects: box_map_batch(fv, rects, padding))
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    return {"boxes": map_graph.num_vertices(),
            "sets": morse_graph.num_vertices()}


def run_data(CMGDB, subdiv, num_points=None, conley=False, batch=False):
    if USE_TORCH:
        raise Skip("torch not applicable to data-defined maps")
    import numpy as np
    lower, upper = LESLIE_DOMAINS[2]
    if num_points is None:
        data = np.loadtxt(EXAMPLES / "Leslie_500.dat")
        X, Y = data[:, :2], data[:, 2:]
    else:
        rng = np.random.default_rng(42)
        X = rng.uniform(lower, upper, size=(num_points, 2))
        Y = numpy_map("leslie", 2)(X)
    F = CMGDB.BoxMapData(X, Y)
    model = CMGDB.Model(subdiv[0], subdiv[1], lower, upper, F)
    if batch:
        if not (hasattr(F, "batch") and hasattr(model, "set_batch_map")):
            raise Skip("no data-map batch")
        model.set_batch_map(F.batch)
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    return {"boxes": map_graph.num_vertices(),
            "sets": morse_graph.num_vertices()}


def run_interval(CMGDB):
    if USE_TORCH:
        raise Skip("torch not applicable to the C++ interval map")
    lower, upper = LESLIE_DOMAINS[2]
    morse_graph = CMGDB.MorseGraphIntvalMap(20, 26, lower, upper,
                                            [19.6, 23.68], "morse_intval.csv")
    return {"boxes": None, "sets": morse_graph.num_vertices()}


def run_expensive(CMGDB, precomputed=False):
    import numpy as np
    lower, upper = LESLIE_DOMAINS[2]
    f_vec = expensive_leslie2d_vec()

    def f_scalar(x):
        return list(np.asarray(f_vec(np.array([x])))[0])

    def F(rect):
        return CMGDB.BoxMap(f_scalar, rect)

    if precomputed:
        if not hasattr(CMGDB, "PrecomputedBoxMap"):
            raise Skip("no PrecomputedBoxMap")
        F_pre = CMGDB.PrecomputedBoxMap(f_vec, lower, upper, 16,
                                        batch_points=16384)
        model = CMGDB.Model(15, 16, lower, upper, F_pre)
        model.set_batch_map(F_pre.batch)
    else:
        model = CMGDB.Model(15, 16, lower, upper, F)
        if not hasattr(model, "set_batch_map"):
            raise Skip("no set_batch_map")
        model.set_batch_map(lambda rects: box_map_batch(f_vec, rects))
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)
    return {"boxes": map_graph.num_vertices(),
            "sets": morse_graph.num_vertices()}


def S(name, fn):
    return (name, fn)


SCENARIOS = dict([
    S("leslie2d_python", lambda C: run_scalar(
        C, "leslie", 2, *LESLIE_DOMAINS[2], (20, 24))),
    S("leslie3d_python", lambda C: run_scalar(
        C, "leslie", 3, *LESLIE_DOMAINS[3], (21, 24))),
    S("leslie4d_python", lambda C: run_scalar(
        C, "leslie", 4, *LESLIE_DOMAINS[4], (22, 26))),
    S("leslie2d_python_deep", lambda C: run_scalar(
        C, "leslie", 2, *LESLIE_DOMAINS[2], (24, 28))),
    S("leslie3d_python_deep", lambda C: run_scalar(
        C, "leslie", 3, *LESLIE_DOMAINS[3], (24, 27))),
    S("leslie2d_python_batch", lambda C: run_scalar(
        C, "leslie", 2, *LESLIE_DOMAINS[2], (20, 24), batch=True)),
    S("leslie3d_python_batch", lambda C: run_scalar(
        C, "leslie", 3, *LESLIE_DOMAINS[3], (21, 24), batch=True)),
    S("leslie2d_python_deep_batch", lambda C: run_scalar(
        C, "leslie", 2, *LESLIE_DOMAINS[2], (24, 28), batch=True)),
    S("leslie2d_interval", lambda C: run_interval(C)),
    S("leslie2d_conley", lambda C: run_scalar(
        C, "leslie", 2, *LESLIE_DOMAINS[2], (20, 22), conley=True)),
    S("leslie3d_conley", lambda C: run_scalar(
        C, "leslie", 3, *LESLIE_DOMAINS[3], (18, 21), conley=True)),
    S("leslie4d_conley", lambda C: run_scalar(
        C, "leslie", 4, *LESLIE_DOMAINS[4], (16, 20), conley=True)),
    S("leslie2d_data", lambda C: run_data(C, (15, 18))),
    S("leslie2d_data_large", lambda C: run_data(C, (14, 16), num_points=50_000)),
    S("leslie2d_data_large_batch", lambda C: run_data(
        C, (14, 16), num_points=50_000, batch=True)),
    S("leslie2d_data_conley", lambda C: run_data(C, (15, 18), conley=True)),
    S("product2d_small", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (6, 10, 4))),
    S("product2d_small_conley", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (6, 10, 4), conley=True)),
    S("product2d_medium", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (10, 14, 8))),
    S("product2d_medium_batch", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (10, 14, 8), batch=True)),
    S("product2d_adaptive", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (14, 18, 8))),
    S("product2d_uniform", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (12, 12, 12))),
    S("product2d_conley", lambda C: run_scalar(
        C, "product", 2, [0.0] * 2, [1.2] * 2, (14, 18, 8), conley=True)),
    S("product3d_uniform", lambda C: run_scalar(
        C, "product", 3, [0.0] * 3, [1.2] * 3, (15, 15, 15))),
    S("product3d_uniform_deep", lambda C: run_scalar(
        C, "product", 3, [0.0] * 3, [1.2] * 3, (18, 18, 18))),
    S("product3d_deep_adaptive", lambda C: run_scalar(
        C, "product", 3, [0.0] * 3, [1.2] * 3, (27, 30, 22))),
    S("product4d_adaptive", lambda C: run_scalar(
        C, "product", 4, [0.0] * 4, [1.2] * 4, (16, 18, 14))),
    S("product4d_deep", lambda C: run_scalar(
        C, "product", 4, [0.0] * 4, [1.2] * 4, (24, 26, 20))),
    S("henon3d", lambda C: run_scalar(
        C, "henon", 3, [-1.5] * 3, [1.5] * 3, (21, 24), padding=True)),
    S("henon3d_conley", lambda C: run_scalar(
        C, "henon", 3, [-1.5] * 3, [1.5] * 3, (12, 15), padding=True,
        conley=True)),
    S("leslie2d_expensive_live", lambda C: run_expensive(C)),
    S("leslie2d_expensive_pre", lambda C: run_expensive(C, precomputed=True)),
])


# ---------------------------------------------------------------------------
# Child mode: run one scenario in the current environment
# ---------------------------------------------------------------------------

def run_child(scenario, out_path):
    import CMGDB
    if USE_TORCH:
        import torch  # noqa: F401 - pay the one-time import outside the timing
    result = {"scenario": scenario}
    try:
        t0 = time.perf_counter()
        info = SCENARIOS[scenario](CMGDB)
        result["seconds"] = time.perf_counter() - t0
        result.update(info)
    except Skip as exc:
        result["skipped"] = str(exc)
    except Exception as exc:  # noqa: BLE001 - report any failure verbatim
        result["error"] = f"{type(exc).__name__}: {exc}"
    Path(out_path).write_text(json.dumps(result))


# ---------------------------------------------------------------------------
# Parent mode: build both versions and orchestrate
# ---------------------------------------------------------------------------

def build_versions(old_rev, workdir):
    old_src = workdir / "old_src"
    if not old_src.exists():
        subprocess.run(["git", "worktree", "add", "--detach",
                        str(old_src), old_rev], cwd=REPO, check=True)
    envs = {}
    for name, src in (("old", old_src), ("new", REPO)):
        venv = workdir / f"venv_{name}"
        python = venv / "bin" / "python"
        stamp = venv / ".built"
        if not stamp.exists():
            print(f"building {name} ({src}) ...", flush=True)
            subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
            log = workdir / f"build_{name}.log"
            with open(log, "w") as fh:
                proc = subprocess.run(
                    [str(venv / "bin" / "pip"), "install", str(src)],
                    stdout=fh, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                sys.exit(f"build of {name} failed; see {log}")
            stamp.touch()
        envs[name] = python
    return envs


def torch_available(python):
    return subprocess.run([str(python), "-c", "import torch"],
                          capture_output=True).returncode == 0


def run_one(python, scenario, tmp, use_torch=False):
    out = Path(tmp) / f"{scenario}.json"
    cmd = [str(python), str(Path(__file__).resolve()),
           "--child", scenario, "--out", str(out)]
    if use_torch:
        cmd.append("--use-torch")
    try:
        proc = subprocess.run(cmd, cwd=tmp, stdout=subprocess.DEVNULL,
                              stderr=subprocess.PIPE, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {"scenario": scenario, "error": f"timeout after {TIMEOUT}s"}
    if proc.returncode != 0 or not out.exists():
        tail = "\n".join(proc.stderr.strip().splitlines()[-4:])
        return {"scenario": scenario, "error": f"exit {proc.returncode}: {tail}"}
    return json.loads(out.read_text())


def fmt_time(seconds):
    if seconds is None:
        return "n/a"
    if seconds < 10:
        return f"{seconds:.2f}s"
    if seconds < 120:
        return f"{seconds:.1f}s"
    return f"{seconds / 60:.1f}m"


def cell(entry):
    if entry is None or "skipped" in entry:
        return "n/a"
    if "error" in entry:
        return "ERROR"
    return fmt_time(entry["seconds"])


def ratio(a, b):
    if (not a or not b or "seconds" not in a or "seconds" not in b):
        return ""
    return f"{a['seconds'] / b['seconds']:.2f}x"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--old-rev", default=DEFAULT_OLD_REV,
                        help=f"git revision to compare against "
                             f"(default: {DEFAULT_OLD_REV}, pre-optimization)")
    parser.add_argument("--scenarios", nargs="*",
                        help="subset of scenarios (default: all)")
    parser.add_argument("--out", default=str(HERE / "version_compare.md"))
    parser.add_argument("--workdir", default=str(REPO / "build" / "version_compare"))
    parser.add_argument("--no-torch", action="store_true",
                        help="skip the torch column even if torch is available")
    parser.add_argument("--list", action="store_true")
    # child-mode internals
    parser.add_argument("--child", help=argparse.SUPPRESS)
    parser.add_argument("--use-torch", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.child:
        global USE_TORCH
        USE_TORCH = args.use_torch
        run_child(args.child, args.out)
        return

    if args.list:
        for name in SCENARIOS:
            print(name)
        return

    selected = args.scenarios or list(SCENARIOS)
    unknown = [s for s in selected if s not in SCENARIOS]
    if unknown:
        sys.exit(f"unknown scenario(s): {', '.join(unknown)}")

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    envs = build_versions(args.old_rev, workdir)

    with_torch = not args.no_torch and torch_available(envs["new"])
    if not args.no_torch and not with_torch:
        print("torch not importable in the new environment; "
              "skipping the torch column (pip install torch into "
              f"{workdir}/venv_new to enable it)")

    flavors = [("old", envs["old"], False), ("new", envs["new"], False)]
    if with_torch:
        flavors.append(("torch", envs["new"], True))

    results = {}
    for scenario in selected:
        results[scenario] = {}
        for name, python, use_torch in flavors:
            print(f"running {scenario} [{name}] ...", flush=True)
            with tempfile.TemporaryDirectory() as tmp:
                entry = run_one(python, scenario, tmp, use_torch)
            results[scenario][name] = entry
            status = ("SKIP" if "skipped" in entry else
                      "ERROR " + entry["error"][:60] if "error" in entry else
                      fmt_time(entry["seconds"]))
            print(f"  -> {status}", flush=True)

    lines = []
    lines.append("# CMGDB version performance comparison\n")
    lines.append(f"Generated {time.strftime('%Y-%m-%d %H:%M')} on "
                 f"{platform.platform()}.\n")
    lines.append(f"**old** = git revision `{args.old_rev}`; "
                 f"**new** = current working tree"
                 + ("; **torch** = new with the map evaluated by torch "
                    "(float64, CPU)" if with_torch else "") + ".\n")
    lines.append("Single repetition per cell, flavors back-to-back per "
                 "scenario for thermal fairness. `n/a` = the mechanism does "
                 "not exist in that version, or (torch column) the scenario "
                 "has no torch-evaluable map: the C++ interval map and the "
                 "data-defined maps, which evaluate no function. In the "
                 "torch column, scalar scenarios run with a torch-backed "
                 "batched map attached; the Conley index phase always "
                 "evaluates the scalar map. Torch runs in float64 so every "
                 "flavor computes identical Morse sets; float32 evaluation "
                 "is roughly twice as fast again but can perturb box "
                 "covers.\n")
    if with_torch:
        lines.append("| Scenario | Old | New | Torch | Old/New | New/Torch |")
        lines.append("|---|---|---|---|---|---|")
    else:
        lines.append("| Scenario | Old | New | Speedup |")
        lines.append("|---|---|---|---|")
    mismatches = []
    for scenario in selected:
        row = results[scenario]
        old, new = row["old"], row["new"]
        if with_torch:
            trc = row.get("torch")
            lines.append(
                f"| {scenario} | {cell(old)} | {cell(new)} | {cell(trc)} "
                f"| {ratio(old, new)} | {ratio(new, trc)} |")
        else:
            lines.append(f"| {scenario} | {cell(old)} | {cell(new)} "
                         f"| {ratio(old, new)} |")
        sets = {name: e.get("sets") for name, e in row.items()
                if e and "seconds" in e}
        if len(set(sets.values())) > 1:
            mismatches.append(f"- `{scenario}`: Morse set counts differ: {sets}")
    if mismatches:
        lines.append("\n## Result mismatches (investigate!)\n")
        lines.extend(mismatches)
    else:
        lines.append("\nAll flavors computed the same number of Morse sets "
                     "on every scenario they could run.\n")
    Path(args.out).write_text("\n".join(lines))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
