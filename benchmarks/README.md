# CMGDB Benchmarks

A timed, correctness-validated benchmark suite for CMGDB. Every scenario
checks its mathematical output (Morse sets, reachability edges, Conley
indices) against a frozen reference before its timings can be trusted —
a performance change that alters the computed dynamics fails loudly.

## Running

```bash
python benchmarks/benchmark.py              # quick suite (~1 min, 3 reps each)
python benchmarks/benchmark.py --heavy      # heavy suite (~3 min, 1 rep each)
python benchmarks/benchmark.py --all        # quick + heavy
python benchmarks/benchmark.py --chafee     # Chafee-Infante 3D suite (opt-in, large)
python benchmarks/benchmark.py --list       # list scenarios
python benchmarks/benchmark.py --scenario leslie2d_python --repeat 5
```

Exit code is nonzero if any scenario fails validation, so the suite can
gate CI.

## Scenarios

The suite covers the main CMGDB usage modes on Leslie population models:

| Scenario | Suite | Mode exercised |
|----------|-------|----------------|
| `leslie2d_python` | quick | Python box map (`BoxMap` corners), 2D |
| `leslie3d_python` | quick | Python box map, 3D |
| `leslie2d_data` | quick | Data-driven `BoxMapData` (`examples/Leslie_500.dat`) |
| `leslie2d_interval` | quick | Built-in C++ interval map (no Python callback) |
| `leslie2d_conley` | quick | Conley-Morse graph (homology computation) |
| `leslie2d_python_deep` | heavy | 2D Python map at deep subdivision |
| `leslie3d_python_deep` | heavy | 3D Python map at deep subdivision |
| `leslie4d_python` | heavy | 4D Python map |
| `leslie2d_data_large` | heavy | `BoxMapData` with 50k synthetic points |

All scenarios are deterministic (corner-mode box maps, seeded data
generation), so references are exact.

## The Chafee-Infante 3D suite

The `chafee` suite (opt-in via `--chafee`; excluded from `--all`) runs the
five computations of a real research workload: the learned 3D
latent-dynamics map of the Chafee-Infante PDE, a tanh MLP
(3 &rarr; 32 &rarr; 32 &rarr; 3) whose trained weights and phase-space
bounds ship in `chafee3d_latent_map.npz` (12 KB). The box map takes
corner images with width padding, evaluated live in vectorized float64
NumPy.

| Scenario | Computation |
|----------|-------------|
| `chafee3d_uniform_16` | Conley-Morse graph, uniform subdiv 16 (65k cells) |
| `chafee3d_uniform_18` | Conley-Morse graph, uniform subdiv 18 (262k cells) |
| `chafee3d_adaptive_18_20_22` | Conley-Morse graph, init 18, subdiv 20/22 |
| `chafee3d_adaptive_21_24_33` | Conley-Morse graph, init 21, subdiv 24/33 (the full-depth adaptive computation) |
| `chafee3d_uniform_24` | Morse graph only, uniform subdiv 24 (16.7M cells, ~1e9 edges, **needs ~10 GB RAM**) |

This is the workload class where the fiber-preboundary and edge-array
reservation behavior matter most (large box images &rarr; large Conley
fibers; a multi-gigabyte CSR edge array at subdiv 24), which the Leslie
and Henon scenarios do not reach at benchmark depths. Total runtime is
roughly 10-15 minutes.

## Reading the results

- **wall (median)** — median wall time over the repetitions.
- **map%** — share of wall time spent inside the Python map function
  (measured by wrapping it). This splits every run into "Python map" vs
  "C++ core" cost, which is the diagnostic for map-evaluation
  optimizations. The interval scenario is 0% by construction.
- **boxes** — final phase-space grid size (where available).
- **MB** — peak RSS of the scenario's subprocess.

Each scenario runs in its own subprocess (clean memory measurement, C++
stdout chatter silenced, one crash cannot kill the suite) from a temporary
working directory (output files like `morse_intval.csv` do not litter the
repo).

## Comparing runs

```bash
python benchmarks/benchmark.py --all --json before.json
# ... apply an optimization, rebuild/reinstall ...
python benchmarks/benchmark.py --all --json after.json
python benchmarks/benchmark.py --compare before.json after.json
```

Two result files are committed in this directory: `baseline_<gitsha>.json`
(timings of the original, unoptimized code) and `optimized_<gitsha>.json`
(timings after the optimization work) — comparing them shows the overall
effect. Machine details are in each file's `meta` block; timings are only
comparable on the same machine. Transient comparison runs named
`after_*.json` are gitignored.

Note that the harness benchmarks the **installed** `CMGDB` package — after
changing C++ code, reinstall (`pip install .`) before benchmarking.

## Comparing against an older revision

`compare_versions.py` builds the current working tree and an older git
revision into two isolated virtual environments, runs an identical
scenario suite under both (feature-detecting mechanisms the old revision
lacks), and writes a Markdown speedup table:

```bash
python benchmarks/compare_versions.py                  # vs 857ec8b (pre-optimization)
python benchmarks/compare_versions.py --old-rev <rev>  # vs any revision
python benchmarks/compare_versions.py --scenarios leslie2d_python henon3d
```

Builds and the old-revision worktree are cached under
`build/version_compare/` (gitignored). Note that old revisions using the
pre-vendoring build require a system sdsl-lite v2 installation to compile.

## Updating references

```bash
python benchmarks/benchmark.py --all --update-refs
```

This overwrites `references.json` for the selected scenarios. It is only
legitimate when the mathematical output is *intended* to change (new
scenario, changed scenario parameters). If an optimization changes a
reference, that optimization has a bug — the entire point of the gate.

`--no-validate` skips the reference check (for calibrating new scenarios);
`--refs PATH` points at an alternative references file.
