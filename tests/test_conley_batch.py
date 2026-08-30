"""Tests for the batched Conley-index phase.

The Conley phase gathers its map evaluations (the index-pair construction
and the top cells of each relative complex) and runs them through the
model's batched evaluator when one is attached, evaluating each unique
rectangle exactly once. These tests pin both the results (identical Morse
graphs and Conley indices with and without a batch map, on uniform and
adaptive runs and through PrecomputedBoxMap) and the behavior (the scalar
map is never called when a batch map is attached, and the scalar fallback
evaluates exactly the set of rectangles the batched path does).

The test map uses only IEEE-exact elementwise operations (division), so the
scalar and vectorized evaluations are bit-identical.
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


class CountingScalar:
    def __init__(self):
        self.calls = 0

    def __call__(self, rect):
        self.calls += 1
        return F_scalar(rect)


class CountingBatch:
    def __init__(self):
        self.calls = 0
        self.rows = 0

    def __call__(self, rects):
        self.calls += 1
        self.rows += np.asarray(rects).shape[0]
        return F_batch(rects)


LOWER, UPPER = [0, 0], [1.2, 1.2]


def build_model(F, subdiv_min=6, subdiv_max=10, subdiv_init=4):
    # subdiv_min < subdiv_max: the final grid is adaptive, so the exit
    # collars of the Morse sets contain coarser cells and the relative
    # complexes refine them (the mixed-depth case).
    return CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, 10000,
                       LOWER, UPPER, F)


def conley_signature(morse_graph, map_graph):
    """Canonical description of the computed dynamics and Conley indices,
    invariant under Morse graph vertex relabeling."""
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
        "conley": [tuple(morse_graph.annotations(v)) for v in order],
    }


@pytest.fixture(scope="module")
def scalar_reference():
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(build_model(F_scalar))
    return conley_signature(morse_graph, map_graph)


def test_scalar_reference_output(scalar_reference):
    # Same dynamics as tests/test_batch_map.py: four Morse sets, and every
    # Conley index defined.
    assert len(scalar_reference["sizes"]) == 4
    assert all("Undefined" not in level
               for annotation in scalar_reference["conley"]
               for level in annotation)


def test_batched_matches_scalar_adaptive(scalar_reference):
    model = build_model(F_scalar)
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(model)
    assert conley_signature(morse_graph, map_graph) == scalar_reference


def test_batched_matches_scalar_small_chunks(scalar_reference):
    # Chunk boundaries inside the Conley phase's gathered evaluations must
    # not affect the result.
    model = build_model(F_scalar)
    model.set_batch_map(F_batch)
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(
        model, batch_chunk_size=100)
    assert conley_signature(morse_graph, map_graph) == scalar_reference


def test_batched_matches_scalar_uniform(scalar_reference):
    del scalar_reference  # uniform run has its own reference
    scalar_model = build_model(F_scalar, subdiv_min=8, subdiv_max=8,
                               subdiv_init=6)
    reference = conley_signature(*CMGDB.ComputeConleyMorseGraph(scalar_model))
    model = build_model(F_scalar, subdiv_min=8, subdiv_max=8, subdiv_init=6)
    model.set_batch_map(F_batch)
    assert conley_signature(*CMGDB.ComputeConleyMorseGraph(model)) == reference


def test_precomputed_box_map_conley(scalar_reference):
    del scalar_reference  # uniform run has its own reference
    scalar_model = build_model(F_scalar, subdiv_min=8, subdiv_max=8,
                               subdiv_init=6)
    reference = conley_signature(*CMGDB.ComputeConleyMorseGraph(scalar_model))
    F_pre = CMGDB.PrecomputedBoxMap(f_vec, LOWER, UPPER, 8)
    model = build_model(F_pre, subdiv_min=8, subdiv_max=8, subdiv_init=6)
    model.set_batch_map(F_pre.batch)
    assert conley_signature(*CMGDB.ComputeConleyMorseGraph(model)) == reference


def test_scalar_map_unused_when_batch_attached():
    # With a batch map attached, no phase -- the Morse-graph construction or
    # the Conley phase -- may fall back to per-rectangle scalar evaluation.
    F = CountingScalar()
    model = build_model(F)
    batch = CountingBatch()
    model.set_batch_map(batch)
    CMGDB.ComputeConleyMorseGraph(model)
    assert F.calls == 0
    assert batch.rows > 0


def test_scalar_fallback_evaluates_same_set_as_batched():
    # The Conley phase gathers its evaluation points once; the batched path
    # and the scalar fallback consume the same gathered vectors, so the
    # number of scalar fallback evaluations must equal the number of rows
    # the batched path processes. Isolate the Conley phase by differencing
    # against Morse-graph-only runs of the same model.
    def batch_rows(compute):
        model = build_model(CountingScalar())
        batch = CountingBatch()
        model.set_batch_map(batch)
        compute(model)
        return batch.rows

    def scalar_calls(compute):
        F = CountingScalar()
        compute(build_model(F))
        return F.calls

    conley_rows = batch_rows(CMGDB.ComputeConleyMorseGraph)
    morse_rows = batch_rows(CMGDB.ComputeMorseGraph)
    conley_calls = scalar_calls(CMGDB.ComputeConleyMorseGraph)
    morse_calls = scalar_calls(CMGDB.ComputeMorseGraph)
    assert conley_rows > morse_rows
    assert conley_rows - morse_rows == conley_calls - morse_calls


def test_conley_index_for_cells_batched():
    model = build_model(F_scalar)
    morse_graph, _ = CMGDB.ComputeConleyMorseGraph(model)
    batched_model = build_model(F_scalar)
    batched_model.set_batch_map(F_batch)
    for v in range(morse_graph.num_vertices()):
        index = CMGDB.ComputeConleyIndexForCells(
            batched_model, morse_graph, morse_graph.morse_set(v))
        assert index == list(morse_graph.annotations(v)), f"node {v}"
