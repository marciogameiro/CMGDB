"""Tests for the cached transition graph and its query API.

Covers the cache_map_graph kwarg and the build_cache() upgrade path, the
automatic edge-reservation kwargs (results must be identical with the
reservation forced on, off, or set explicitly), the native cached-graph
reachability queries (validated against brute-force BFS oracles over the
returned map_graph), and ComputeConleyIndexForCells (validated against the
per-Morse-set Conley annotations of ComputeConleyMorseGraph).
"""

import math

import pytest

import CMGDB


def product_model(dim=2, subdiv=10):
    """Uniform-grid product map x -> x/(2-x): 5^dim Morse sets, known DAG."""
    def f(x):
        return [x[d] / (2.0 - x[d]) for d in range(dim)]

    def F(rect):
        return CMGDB.BoxMap(f, rect)

    return CMGDB.Model(subdiv, subdiv, subdiv, 10000,
                       [0.0] * dim, [1.2] * dim, F)


def leslie_model(subdiv_min=12, subdiv_max=14):
    """Adaptive 2D Leslie model: the final grid has cells at mixed depths."""
    theta = [19.6, 23.68]

    def f(x):
        s = x[0] + x[1]
        return [(theta[0] * x[0] + theta[1] * x[1]) * math.exp(-0.1 * s),
                0.7 * x[0]]

    def F(rect):
        return CMGDB.BoxMap(f, rect)

    return CMGDB.Model(subdiv_min, subdiv_max,
                       [-0.001, -0.001], [90.0, 70.0], F)


def all_adjacencies(map_graph):
    return [list(map_graph.adjacencies(v))
            for v in range(map_graph.num_vertices())]


# ---------------------------------------------------------------------------
# cache_map_graph / build_cache
# ---------------------------------------------------------------------------

def test_returned_map_graph_is_lazy_by_default():
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(product_model())
    assert not map_graph.has_cache()
    assert map_graph.num_cached_edges() == 0
    # Lazy adjacency queries still work.
    assert map_graph.num_vertices() > 0
    assert isinstance(map_graph.adjacencies(0), list)


def test_cache_map_graph_returns_cached_graph_with_same_adjacencies():
    mg_lazy, gr_lazy = CMGDB.ComputeMorseGraph(product_model())
    mg_eager, gr_eager = CMGDB.ComputeMorseGraph(product_model(),
                                                 cache_map_graph=True)
    assert gr_eager.has_cache()
    assert gr_eager.num_cached_edges() > 0
    assert mg_lazy.num_vertices() == mg_eager.num_vertices()
    assert all_adjacencies(gr_lazy) == all_adjacencies(gr_eager)


def test_build_cache_upgrades_lazy_graph_and_is_idempotent():
    _, map_graph = CMGDB.ComputeMorseGraph(product_model())
    lazy = all_adjacencies(map_graph)
    map_graph.build_cache()
    assert map_graph.has_cache()
    edges = map_graph.num_cached_edges()
    assert edges > 0
    map_graph.build_cache()  # no-op
    assert map_graph.num_cached_edges() == edges
    assert all_adjacencies(map_graph) == lazy


# ---------------------------------------------------------------------------
# Edge reservation kwargs
# ---------------------------------------------------------------------------

def reference_of(morse_graph):
    sets = [sorted(morse_graph.morse_set(v))
            for v in range(morse_graph.num_vertices())]
    return sorted(sets), sorted(morse_graph.edges_unreduced())


def test_reserve_kwargs_change_nothing_mathematical():
    baseline_mg, baseline_gr = CMGDB.ComputeMorseGraph(
        product_model(), cache_map_graph=True)
    forced_mg, forced_gr = CMGDB.ComputeMorseGraph(
        product_model(), cache_map_graph=True, reserve_min_edges=1)
    explicit_mg, explicit_gr = CMGDB.ComputeMorseGraph(
        product_model(), cache_map_graph=True,
        reserve_edges=100_000, reserve_min_edges=1)
    assert reference_of(baseline_mg) == reference_of(forced_mg)
    assert reference_of(baseline_mg) == reference_of(explicit_mg)
    assert (all_adjacencies(baseline_gr) == all_adjacencies(forced_gr)
            == all_adjacencies(explicit_gr))
    assert (baseline_gr.num_cached_edges() == forced_gr.num_cached_edges()
            == explicit_gr.num_cached_edges())


# ---------------------------------------------------------------------------
# Native cached-graph reachability queries vs brute-force BFS oracles
# ---------------------------------------------------------------------------

def bfs_reach(adjacency, seeds):
    seen = set(seeds)
    stack = list(seeds)
    while stack:
        v = stack.pop()
        for w in adjacency[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return seen


@pytest.fixture(scope="module")
def cached_run():
    model = product_model()
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model,
                                                     cache_map_graph=True)
    adjacency = all_adjacencies(map_graph)
    morse_sets = [set(morse_graph.morse_set(v))
                  for v in range(morse_graph.num_vertices())]
    return model, morse_graph, map_graph, adjacency, morse_sets


def test_reachability_masks_match_bfs_oracle(cached_run):
    _, morse_graph, map_graph, adjacency, morse_sets = cached_run
    n = map_graph.num_vertices()
    queries = list(range(n))
    masks = CMGDB.MorseReachabilityMasks(map_graph, morse_graph, queries)
    assert len(masks) == n
    for q in queries:
        reach = bfs_reach(adjacency, [q])
        expected = 0
        for node, cells in enumerate(morse_sets):
            if reach & cells:
                expected |= 1 << node
        assert masks[q] == expected, f"mask mismatch at cell {q}"


def test_singleton_reachability_matches_bfs_oracle(cached_run):
    _, morse_graph, map_graph, adjacency, morse_sets = cached_run
    n = map_graph.num_vertices()
    queries = list(range(n))
    summary = CMGDB.MorseSingletonReachability(map_graph, morse_graph, queries)
    for q in queries:
        reach = bfs_reach(adjacency, [q])
        nodes = {node for node, cells in enumerate(morse_sets)
                 if reach & cells}
        if not nodes:
            expected = -1
        elif len(nodes) == 1:
            expected = next(iter(nodes))
        else:
            expected = -2
        assert summary[q] == expected, f"summary mismatch at cell {q}"


def test_directed_path_cells_match_bfs_oracle(cached_run):
    _, morse_graph, map_graph, adjacency, morse_sets = cached_run
    n = map_graph.num_vertices()
    reverse = [[] for _ in range(n)]
    for v, targets in enumerate(adjacency):
        for w in targets:
            reverse[w].append(v)
    # Pick a source node with outgoing Morse edges and a target it reaches.
    edges = [e for e in morse_graph.edges_unreduced() if e[0] != e[1]]
    assert edges, "product map should have nontrivial Morse edges"
    source_node, target_node = edges[0]
    result = CMGDB.MorseDirectedPathCells(
        map_graph, morse_graph, [source_node], [target_node])
    forward = bfs_reach(adjacency, sorted(morse_sets[source_node]))
    backward = bfs_reach(reverse, sorted(morse_sets[target_node]))
    expected = sorted(forward & backward)
    assert sorted(result) == expected


def test_native_queries_require_cached_graph():
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(product_model())
    assert not map_graph.has_cache()
    with pytest.raises(RuntimeError, match="cached MapGraph"):
        CMGDB.MorseReachabilityMasks(map_graph, morse_graph, [0])
    with pytest.raises(RuntimeError, match="cached MapGraph"):
        CMGDB.MorseSingletonReachability(map_graph, morse_graph, [0])
    with pytest.raises(RuntimeError, match="cached MapGraph"):
        CMGDB.MorseDirectedPathCells(map_graph, morse_graph, [0], [0])


# ---------------------------------------------------------------------------
# ComputeConleyIndexForCells
# ---------------------------------------------------------------------------

def test_conley_index_for_cells_matches_annotations_uniform():
    model = product_model(subdiv=8)
    morse_graph, _ = CMGDB.ComputeConleyMorseGraph(model)
    for v in range(morse_graph.num_vertices()):
        index = CMGDB.ComputeConleyIndexForCells(
            model, morse_graph, morse_graph.morse_set(v))
        assert index == list(morse_graph.annotations(v)), f"node {v}"


def test_conley_index_for_cells_matches_annotations_adaptive():
    # Adaptive run: the final joined grid holds cells at mixed depths, so
    # this exercises the common-refinement path of the chomp machinery.
    model = leslie_model()
    morse_graph, _ = CMGDB.ComputeConleyMorseGraph(model)
    assert morse_graph.num_vertices() > 0
    for v in range(morse_graph.num_vertices()):
        index = CMGDB.ComputeConleyIndexForCells(
            model, morse_graph, morse_graph.morse_set(v))
        assert index == list(morse_graph.annotations(v)), f"node {v}"


def test_conley_index_for_cells_validates_input():
    model = product_model(subdiv=8)
    morse_graph, _ = CMGDB.ComputeConleyMorseGraph(model)
    with pytest.raises(Exception):
        CMGDB.ComputeConleyIndexForCells(model, morse_graph, [2**62])
