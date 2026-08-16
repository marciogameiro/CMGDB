# Vendored dependency: sdsl-lite v3.0.4 (header-only)

Source: https://github.com/xxsds/sdsl-lite, tag `v3.0.4`
(release tarball SHA256 `9bade92986d5b6dae15b723f6b2d87b14842e56295558f88c8daaeb33c46967e`).

This is the community-maintained, header-only, BSD-3-Clause fork of the
Succinct Data Structure Library (the original `simongog/sdsl-lite` v2 is a
GPLv3 compiled library). Only the `include/` tree plus `LICENSE` and
`AUTHORS` are vendored; CMGDB uses it from `SuccinctTree.h` and
`RankSelect.h`, which back the optional succinct grid and are compiled
only when `CMGDB_USE_SUCCINCT` is defined (default builds use PointerGrid
and never include sdsl). Note: these headers use GCC/Clang extensions
(`__attribute__`, `cxxabi.h`) and do not compile under MSVC.

Vendored rather than fetched or installed so that source builds
(`pip install CMGDB`) are self-contained, need no network at build time,
and cannot pick up an incompatible sdsl v2 from the system include paths.

To upgrade: replace `include/` with the new release's `include/`, update the
tag and checksum here, rebuild, and run `pytest tests/` plus
`python benchmarks/benchmark.py --all` (the frozen references validate that
the succinct grid structures still behave identically).
