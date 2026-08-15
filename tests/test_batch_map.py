"""Equivalence tests for the transition-graph cache and batched map evaluation.

The test map uses only IEEE-exact elementwise operations (division), so the
scalar and vectorized evaluations are bit-identical and every computation
below must produce exactly the same Morse graph.
"""

import numpy as np
import pytest
import CMGDB


def f_vec(X):
    return np.column_stack([X[:, 0] / (2.0 - X[:, 0]), X[:, 1] / (2.0 - X[:, 1])])


def f_scalar(x):
    return list(f_vec(np.array([x]))[0])


def F_scalar(rect):
    return CMGDB.BoxMap(f_scalar, rect, padding=False)


def F_batch(rects):
    return CMGDB.BoxMapBatch(f_vec, rects, padding=False)


def build_model():
    subdiv_min, subdiv_max = 6, 10
    subdiv_init, subdiv_limit = 4, 10000
    lower_bounds, upper_bounds = [0, 0], [1.2, 1.2]
    return CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, subdiv_limit,
                       lower_bounds, upper_bounds, F_scalar)


def signature(morse_graph, map_graph):
    """Canonical description of the computed dynamics, invariant under
    Morse graph vertex relabeling."""
    num_vertices = morse_graph.num_vertices()
    morse_sets = [sorted(morse_graph.morse_set(v)) for v in range(num_vertices)]
    order = sorted(range(num_vertices), key=lambda v: morse_sets[v][0])
    relabel = {v: i for i, v in enumerate(order)}
    return {
        "sizes": [len(morse_sets[v]) for v in order],
        "min_box": [morse_sets[v][0] for v in order],
        "edges": sorted((relabel[u], relabel[w])
                        for u, w in morse_graph.edges_unreduced()),
        "phase_size": map_graph.num_vertices(),
    }


@pytest.fixture(scope="module")
def reference_signature():
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(build_model())
    return signature(morse_graph, map_graph)


def test_reference_output(reference_signature):
    # Same computation as tests/test_basic.py: four Morse sets
    assert len(reference_signature["sizes"]) == 4


def test_batch_map_matches_scalar(reference_signature):
    model = build_model()
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)
    assert signature(morse_graph, map_graph) == reference_signature


def test_cache_disabled_matches(reference_signature):
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(
        build_model(), cache_transition_graph=False)
    assert signature(morse_graph, map_graph) == reference_signature


def test_small_chunks_match(reference_signature):
    # Chunk boundaries must not affect the result
    model = build_model()
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model, batch_chunk_size=100)
    assert signature(morse_graph, map_graph) == reference_signature


def test_max_cached_edges_fallback_matches(reference_signature):
    # A tiny edge cap forces the cache to be abandoned at every level;
    # the on-demand fallback must produce the identical result
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(
        build_model(), max_cached_edges=10)
    assert signature(morse_graph, map_graph) == reference_signature


def test_batch_map_with_conley_index(reference_signature):
    model = build_model()
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(model)
    assert signature(morse_graph, map_graph) == reference_signature
    annotations = sorted(tuple(morse_graph.annotations(v))
                         for v in range(morse_graph.num_vertices()))
    assert annotations == [('0', '0', 'x-1'), ('0', 'x-1', '0'),
                           ('0', 'x-1', '0'), ('x-1', '0', '0')]


def test_wrong_batch_shape_raises():
    model = build_model()
    model.set_batch_map(lambda rects: np.zeros((3, 3)))
    with pytest.raises(RuntimeError):
        CMGDB.ComputeMorseGraph(model)


def test_non_array_batch_return_raises():
    model = build_model()
    model.set_batch_map(lambda rects: "not an array")
    with pytest.raises(RuntimeError):
        CMGDB.ComputeMorseGraph(model)


def test_set_batch_map_requires_function_map():
    # A model built without a map function F cannot take a batch map
    model = CMGDB.Model(6, 10, [0, 0], [1.2, 1.2])
    with pytest.raises(Exception):
        model.set_batch_map(F_batch)


def test_boxmap_batch_matches_boxmap_rowwise():
    rng = np.random.default_rng(7)
    lows = rng.uniform(0.0, 0.5, size=(50, 2))
    ups = lows + rng.uniform(0.01, 0.5, size=(50, 2))
    rects = np.hstack([lows, ups])
    for padding in (False, True):
        batched = CMGDB.BoxMapBatch(f_vec, rects, padding=padding)
        for i in range(rects.shape[0]):
            scalar = CMGDB.BoxMap(f_scalar, list(rects[i]), padding=padding)
            assert np.array_equal(batched[i], np.array(scalar))
