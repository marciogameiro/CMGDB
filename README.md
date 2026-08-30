# CMGDB
Conley Morse Graph Database

## Overview

This project uses combinatorial and topological methods to compute dynamics of discrete dynamical systems.

## Installation

Install the latest tagged version:

	pip install CMGDB

To uninstall:

	pip uninstall CMGDB

## Documentation and examples

To get started on how to run the code see the examples in the Jupyter notebooks in the [examples](examples) folder.

In particular the notebooks [Examples.ipynb](examples/Examples.ipynb), [Gaussian\_Process\_Example.ipynb](examples/Gaussian_Process_Example.ipynb), and [Conley\_Index\_Examples.ipynb](examples/Conley_Index_Examples.ipynb) present basic examples on how to run the code and are a good starting point.

Here is an old [survey](http://chomp.rutgers.edu/Projects/survey/cmdbSurvey.pdf) and a
[talk](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf) that might be useful.

## Performance options and the transition-graph cache

`ComputeMorseGraph` and `ComputeConleyMorseGraph` return a pair `(morse_graph, map_graph)` and accept keyword arguments controlling the transition-graph machinery (the defaults are right for most runs):

* `cache_transition_graph` (default `True`) — cache the per-level transition graph used internally by the SCC/reachability passes, halving the map evaluations per subdivision level. Set `False` for a memory-lean run that re-evaluates the map on demand.
* `batch_chunk_size` (default `65536`) — rectangles per batched map call when a batch map is attached with `model.set_batch_map` (`0` means one call for the whole grid). The Conley-index phase of `ComputeConleyMorseGraph` (and `ComputeConleyIndexForCells`) also gathers its map evaluations into these chunks, evaluating each rectangle exactly once — without a batch map the evaluations are scalar but still deduplicated, so attaching a batch map speeds up every phase, not just the transition graph.
* `max_cached_edges` (default `0` = unlimited) — abandon a cache measured (or confidently projected) to exceed this many edges and fall back to on-demand evaluation.
* `reserve_edges` / `reserve_min_edges` (defaults `0` / `2**24`) — up-front sizing of the flat edge array. By default the final edge count is projected from the first chunk and twice that is reserved, which avoids the reallocation spikes of multi-gigabyte graphs on deep grids; a positive `reserve_edges` reserves exactly that many instead. Reservation only engages once the projection reaches `reserve_min_edges`.
* `cache_map_graph` (default `False`) — eagerly cache the **returned** `map_graph` (one extra full batched map pass over the final grid, after which its adjacency queries are O(1) array lookups). The default returns a lazy `map_graph` that evaluates the map per `adjacencies` query; `map_graph.build_cache()` upgrades it to the cached form later. `map_graph.has_cache()` and `map_graph.num_cached_edges()` report the state.

A cached `map_graph` also unlocks the native post-processing queries (all of which release the GIL and refuse a lazy graph):

* `MorseReachabilityMasks(map_graph, morse_graph, cells)` — for each queried cell, a `uint64` bitmask of the Morse nodes reachable through the box dynamics (bit `i` = Morse node `i`); the basin-of-attraction primitive.
* `MorseSingletonReachability(map_graph, morse_graph, cells)` — per cell, the single reachable Morse node id, `-1` for none, `-2` for several.
* `MorseDirectedPathCells(map_graph, morse_graph, sources, targets)` — the cells lying on some directed path from the source Morse nodes to the target Morse nodes (candidate connecting-orbit regions).
* `ComputeConleyIndexForCells(model, morse_graph, cells)` — the homological Conley index of an arbitrary cell subset of the final grid (cells of different sizes from adaptive runs are handled exactly, by common refinement).

The C++ core is silent by default; rebuild with the `CMG_VERBOSE` preprocessor define (uncomment it at the top of `src/CMGDB/_cmgdb/CMGDB.cpp`) to restore the progress and diagnostic prints.

## Installing from source and dependencies

To install from source you need a C++ compiler and [Boost](https://www.boost.org/) installed. Assuming you have Boost installed in your system, you can install from source with the command:

	pip install --force-reinstall --no-deps --no-cache-dir git+https://github.com/marciogameiro/CMGDB.git

Alternatively, you can clone the GitHub repository and install with:

	git clone https://github.com/marciogameiro/CMGDB.git
	cd CMGDB
	pip install .
