#!/usr/bin/env python3
"""Compare CMGDB performance between two versions of the code.

Builds the current working tree and an older git revision into two isolated
virtual environments, runs an identical scenario suite under both, and
writes a Markdown speedup table. Mechanisms that do not exist in the old
revision (batched maps, data-map batch, precomputed box maps) are
feature-detected and reported as n/a there.

Usage:
  python benchmarks/compare_versions.py                  # vs 857ec8b (pre-optimization)
  python benchmarks/compare_versions.py --old-rev v1.3.2 # vs any git revision
  python benchmarks/compare_versions.py --scenarios leslie2d_python henon3d
  python benchmarks/compare_versions.py --out my_table.md

The two builds take a few minutes; results land in version_compare.md (or
--out) next to this script. Scenarios run one repetition each, old and new
back-to-back per scenario for thermal fairness. The scenario suite mirrors
benchmarks/benchmark.py but validates only that both versions compute the
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

# ---------------------------------------------------------------------------
# Maps (defined here so both versions run identical code)
# ---------------------------------------------------------------------------

LESLIE_THETA = [19.6, 23.68, 19.6, 23.68]
LESLIE_DOMAINS = {
    2: ([-0.001, -0.001], [90.0, 70.0]),
    3: ([-0.001, -0.001, -0.001], [90.0, 65.0, 45.0]),
    4: ([-0.001, -0.001, -0.001, -0.001], [90.0, 65.0, 45.0, 32.0]),
}


def leslie_scalar(dim):
    theta = LESLIE_THETA[:dim]

    def f(x):
        s = sum(x)
        y0 = sum(t * xi for t, xi in zip(theta, x)) * math.exp(-0.1 * s)
        return [y0] + [0.7 * xi for xi in x[:dim - 1]]

    return f


def leslie_vec(dim):
    import numpy as np
    theta = np.array(LESLIE_THETA[:dim])

    def f(X):
        s = X.sum(axis=1)
        y0 = (X @ theta) * np.exp(-0.1 * s)
        return np.column_stack([y0] + [0.7 * X[:, i] for i in range(dim - 1)])

    return f


def product_scalar(dim):
    def f(x):
        return [x[d] / (2.0 - x[d]) for d in range(dim)]

    return f


def product_vec(dim):
    def f(X):
        return X / (2.0 - X)

    return f


def henon3d_scalar(x):
    return [1.0 - 1.4 * x[0] * x[0] + 0.3 * x[1], x[0], x[1]]


def expensive_leslie2d_vec(width=1024, seed=7, amplitude=0.5):
    import numpy as np
    rng = np.random.default_rng(seed)
    W1 = rng.normal(0.0, 1.0, (2, width))
    b1 = rng.normal(0.0, 1.0, width)
    W2 = rng.normal(0.0, 1.0 / np.sqrt(width), (width, 2))
    b2 = rng.normal(0.0, 1.0, 2)
    scale_in = np.array([45.0, 35.0])
    base = leslie_vec(2)

    def f(P):
        hidden = np.tanh(P / scale_in @ W1 + b1)
        return base(P) + amplitude * np.tanh(hidden @ W2 + b2)

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


def run_scalar(CMGDB, scalar_f, lower, upper, subdiv, conley=False,
               padding=False, batch_vec=None):
    def F(rect):
        return CMGDB.BoxMap(scalar_f, rect, padding=padding)

    if len(subdiv) == 2:
        model = CMGDB.Model(subdiv[0], subdiv[1], lower, upper, F)
    else:
        model = CMGDB.Model(subdiv[0], subdiv[1], subdiv[2], 10000,
                            lower, upper, F)
    if batch_vec is not None:
        if not hasattr(model, "set_batch_map"):
            raise Skip("no set_batch_map")
        model.set_batch_map(lambda rects: box_map_batch(batch_vec, rects, padding))
    compute = CMGDB.ComputeConleyMorseGraph if conley else CMGDB.ComputeMorseGraph
    morse_graph, map_graph = compute(model)
    return {"boxes": map_graph.num_vertices(),
            "sets": morse_graph.num_vertices()}


def run_data(CMGDB, subdiv, num_points=None, conley=False, batch=False):
    import numpy as np
    lower, upper = LESLIE_DOMAINS[2]
    if num_points is None:
        data = np.loadtxt(EXAMPLES / "Leslie_500.dat")
        X, Y = data[:, :2], data[:, 2:]
    else:
        rng = np.random.default_rng(42)
        X = rng.uniform(lower, upper, size=(num_points, 2))
        Y = leslie_vec(2)(X)
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
    lower, upper = LESLIE_DOMAINS[2]
    morse_graph = CMGDB.MorseGraphIntvalMap(20, 26, lower, upper,
                                            [19.6, 23.68], "morse_intval.csv")
    return {"boxes": None, "sets": morse_graph.num_vertices()}


def run_expensive(CMGDB, precomputed=False):
    import numpy as np
    lower, upper = LESLIE_DOMAINS[2]
    f_vec = expensive_leslie2d_vec()

    def f_scalar(x):
        return list(f_vec(np.array([x]))[0])

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
        C, leslie_scalar(2), *LESLIE_DOMAINS[2], (20, 24))),
    S("leslie3d_python", lambda C: run_scalar(
        C, leslie_scalar(3), *LESLIE_DOMAINS[3], (21, 24))),
    S("leslie4d_python", lambda C: run_scalar(
        C, leslie_scalar(4), *LESLIE_DOMAINS[4], (22, 26))),
    S("leslie2d_python_deep", lambda C: run_scalar(
        C, leslie_scalar(2), *LESLIE_DOMAINS[2], (24, 28))),
    S("leslie3d_python_deep", lambda C: run_scalar(
        C, leslie_scalar(3), *LESLIE_DOMAINS[3], (24, 27))),
    S("leslie2d_python_batch", lambda C: run_scalar(
        C, leslie_scalar(2), *LESLIE_DOMAINS[2], (20, 24),
        batch_vec=leslie_vec(2))),
    S("leslie3d_python_batch", lambda C: run_scalar(
        C, leslie_scalar(3), *LESLIE_DOMAINS[3], (21, 24),
        batch_vec=leslie_vec(3))),
    S("leslie2d_python_deep_batch", lambda C: run_scalar(
        C, leslie_scalar(2), *LESLIE_DOMAINS[2], (24, 28),
        batch_vec=leslie_vec(2))),
    S("leslie2d_interval", lambda C: run_interval(C)),
    S("leslie2d_conley", lambda C: run_scalar(
        C, leslie_scalar(2), *LESLIE_DOMAINS[2], (20, 22), conley=True)),
    S("leslie3d_conley", lambda C: run_scalar(
        C, leslie_scalar(3), *LESLIE_DOMAINS[3], (18, 21), conley=True)),
    S("leslie4d_conley", lambda C: run_scalar(
        C, leslie_scalar(4), *LESLIE_DOMAINS[4], (16, 20), conley=True)),
    S("leslie2d_data", lambda C: run_data(C, (15, 18))),
    S("leslie2d_data_large", lambda C: run_data(C, (14, 16), num_points=50_000)),
    S("leslie2d_data_large_batch", lambda C: run_data(
        C, (14, 16), num_points=50_000, batch=True)),
    S("leslie2d_data_conley", lambda C: run_data(C, (15, 18), conley=True)),
    S("product2d_small", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (6, 10, 4))),
    S("product2d_small_conley", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (6, 10, 4), conley=True)),
    S("product2d_medium", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (10, 14, 8))),
    S("product2d_medium_batch", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (10, 14, 8),
        batch_vec=product_vec(2))),
    S("product2d_adaptive", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (14, 18, 8))),
    S("product2d_uniform", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (12, 12, 12))),
    S("product2d_conley", lambda C: run_scalar(
        C, product_scalar(2), [0.0] * 2, [1.2] * 2, (14, 18, 8), conley=True)),
    S("product3d_uniform", lambda C: run_scalar(
        C, product_scalar(3), [0.0] * 3, [1.2] * 3, (15, 15, 15))),
    S("product3d_uniform_deep", lambda C: run_scalar(
        C, product_scalar(3), [0.0] * 3, [1.2] * 3, (18, 18, 18))),
    S("product3d_deep_adaptive", lambda C: run_scalar(
        C, product_scalar(3), [0.0] * 3, [1.2] * 3, (27, 30, 22))),
    S("product4d_adaptive", lambda C: run_scalar(
        C, product_scalar(4), [0.0] * 4, [1.2] * 4, (16, 18, 14))),
    S("product4d_deep", lambda C: run_scalar(
        C, product_scalar(4), [0.0] * 4, [1.2] * 4, (24, 26, 20))),
    S("henon3d", lambda C: run_scalar(
        C, henon3d_scalar, [-1.5] * 3, [1.5] * 3, (21, 24), padding=True)),
    S("henon3d_conley", lambda C: run_scalar(
        C, henon3d_scalar, [-1.5] * 3, [1.5] * 3, (12, 15), padding=True,
        conley=True)),
    S("leslie2d_expensive_live", lambda C: run_expensive(C)),
    S("leslie2d_expensive_pre", lambda C: run_expensive(C, precomputed=True)),
])


# ---------------------------------------------------------------------------
# Child mode: run one scenario in the current environment
# ---------------------------------------------------------------------------

def run_child(scenario, out_path):
    import CMGDB
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


def run_one(python, scenario, tmp):
    out = Path(tmp) / f"{scenario}.json"
    cmd = [str(python), str(Path(__file__).resolve()),
           "--child", scenario, "--out", str(out)]
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--old-rev", default=DEFAULT_OLD_REV,
                        help=f"git revision to compare against "
                             f"(default: {DEFAULT_OLD_REV}, pre-optimization)")
    parser.add_argument("--scenarios", nargs="*",
                        help="subset of scenarios (default: all)")
    parser.add_argument("--out", default=str(HERE / "version_compare.md"))
    parser.add_argument("--workdir", default=str(REPO / "build" / "version_compare"))
    parser.add_argument("--list", action="store_true")
    # child-mode internals
    parser.add_argument("--child", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.child:
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

    results = {}
    for scenario in selected:
        results[scenario] = {}
        for name in ("old", "new"):
            print(f"running {scenario} [{name}] ...", flush=True)
            with tempfile.TemporaryDirectory() as tmp:
                entry = run_one(envs[name], scenario, tmp)
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
                 f"**new** = current working tree.\n")
    lines.append("Single repetition per cell, old and new back-to-back per "
                 "scenario for thermal fairness. `n/a` = the mechanism does "
                 "not exist in that version.\n")
    lines.append("| Scenario | Old | New | Speedup |")
    lines.append("|---|---|---|---|")
    mismatches = []
    for scenario in selected:
        old, new = results[scenario]["old"], results[scenario]["new"]
        speedup = ""
        if "seconds" in old and "seconds" in new:
            speedup = f"{old['seconds'] / new['seconds']:.2f}x"
        lines.append(f"| {scenario} | {cell(old)} | {cell(new)} | {speedup} |")
        if ("seconds" in old and "seconds" in new
                and old.get("sets") != new.get("sets")):
            mismatches.append(
                f"- `{scenario}`: Morse set counts differ: "
                f"old={old.get('sets')} new={new.get('sets')}")
    if mismatches:
        lines.append("\n## Result mismatches (investigate!)\n")
        lines.extend(mismatches)
    else:
        lines.append("\nBoth versions computed the same number of Morse sets "
                     "on every scenario they could run.\n")
    Path(args.out).write_text("\n".join(lines))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
