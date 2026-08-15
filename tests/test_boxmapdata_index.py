"""Equivalence tests: the indexed BoxMapData must reproduce the linear-scan
reference implementation BoxMapDataLinear exactly, on every code path."""

from pathlib import Path

import numpy as np
import pytest
import CMGDB

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"

THETA = np.array([19.6, 23.68])
LOWER_BOUNDS = [-0.001, -0.001]
UPPER_BOUNDS = [90.0, 70.0]


def leslie(X):
    s = X.sum(axis=1)
    y0 = (X @ THETA) * np.exp(-0.1 * s)
    return np.column_stack([y0, 0.7 * X[:, 0]])


@pytest.fixture(scope="module")
def data():
    rng = np.random.default_rng(3)
    X = rng.uniform(LOWER_BOUNDS, UPPER_BOUNDS, size=(2000, 2))
    return X, leslie(X)


def random_rects(rng, count, scale):
    lows = rng.uniform(-10.0, 100.0, size=(count, 2))
    widths = rng.uniform(0.0, scale, size=(count, 2))
    return np.hstack([lows, lows + widths])


def test_map_points_identical_on_random_rects(data):
    X, Y = data
    indexed = CMGDB.BoxMapData(X, Y)
    linear = CMGDB.BoxMapDataLinear(X, Y)
    rng = np.random.default_rng(11)
    # Tiny, medium, huge, and mostly-outside rectangles
    for scale in (0.5, 5.0, 50.0, 500.0):
        for rect in random_rects(rng, 50, scale):
            assert np.array_equal(indexed.map_points(list(rect)),
                                  linear.map_points(list(rect)))


def test_map_points_boundary_points(data):
    # Rectangles whose bounds coincide exactly with data point coordinates:
    # the closed containment test must keep behaving identically
    X, Y = data
    indexed = CMGDB.BoxMapData(X, Y)
    linear = CMGDB.BoxMapDataLinear(X, Y)
    for i in (0, 7, 1999):
        point = X[i]
        degenerate = list(point) + list(point)
        assert np.array_equal(indexed.map_points(degenerate),
                              linear.map_points(degenerate))
        corner = list(point) + list(point + 3.0)
        assert np.array_equal(indexed.map_points(corner),
                              linear.map_points(corner))


def test_compute_identical_across_options(data):
    X, Y = data
    rng = np.random.default_rng(13)
    rects = random_rects(rng, 40, 5.0)
    option_sets = [
        dict(map_empty='interp'),
        dict(map_empty='outside', lower_bounds=LOWER_BOUNDS, upper_bounds=UPPER_BOUNDS),
        dict(map_empty='interp', domain_padding=False),
        dict(map_empty='interp', padding=True),
        dict(map_empty='interp', domain_padding=False, padding=True),
    ]
    for options in option_sets:
        indexed = CMGDB.BoxMapData(X, Y, **options)
        linear = CMGDB.BoxMapDataLinear(X, Y, **options)
        for rect in rects:
            assert indexed.compute(list(rect)) == linear.compute(list(rect))


def test_terminate_mode_raises_identically(data):
    X, Y = data
    empty_rect = [200.0, 200.0, 201.0, 201.0]
    indexed = CMGDB.BoxMapData(X, Y, map_empty='terminate')
    linear = CMGDB.BoxMapDataLinear(X, Y, map_empty='terminate')
    with pytest.raises(ValueError):
        indexed.compute(empty_rect)
    with pytest.raises(ValueError):
        linear.compute(empty_rect)


def test_forced_index_on_small_dataset_identical(data):
    # Below the auto threshold the index is skipped; forcing it on must
    # still give identical results
    X, Y = data
    X_small, Y_small = X[:200], Y[:200]
    forced = CMGDB.BoxMapData(X_small, Y_small, use_index=True)
    assert forced._use_index
    auto = CMGDB.BoxMapData(X_small, Y_small)
    assert not auto._use_index
    linear = CMGDB.BoxMapDataLinear(X_small, Y_small)
    rng = np.random.default_rng(19)
    for rect in random_rects(rng, 50, 10.0):
        expected = linear.map_points(list(rect))
        assert np.array_equal(forced.map_points(list(rect)), expected)
        assert np.array_equal(auto.map_points(list(rect)), expected)


def test_batch_matches_scalar(data):
    X, Y = data
    indexed = CMGDB.BoxMapData(X, Y)
    rng = np.random.default_rng(17)
    rects = random_rects(rng, 30, 5.0)
    assert indexed.batch(rects) == [indexed(list(rect)) for rect in rects]


def morse_signature(morse_graph, map_graph):
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


def test_morse_graph_identical_on_leslie_data():
    dataset = np.loadtxt(EXAMPLES_DIR / "Leslie_500.dat")
    X, Y = dataset[:, :2], dataset[:, 2:]
    signatures = []
    for cls in (CMGDB.BoxMapDataLinear, CMGDB.BoxMapData):
        F = cls(X, Y)
        model = CMGDB.Model(10, 12, LOWER_BOUNDS, UPPER_BOUNDS, F)
        morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)
        signatures.append(morse_signature(morse_graph, map_graph))
    assert signatures[0] == signatures[1]


def test_morse_graph_identical_with_batch():
    dataset = np.loadtxt(EXAMPLES_DIR / "Leslie_500.dat")
    X, Y = dataset[:, :2], dataset[:, 2:]
    F = CMGDB.BoxMapData(X, Y)
    model = CMGDB.Model(10, 12, LOWER_BOUNDS, UPPER_BOUNDS, F)
    reference_mg, reference_g = CMGDB.ComputeMorseGraph(model)
    F_batch = CMGDB.BoxMapData(X, Y)
    model_batch = CMGDB.Model(10, 12, LOWER_BOUNDS, UPPER_BOUNDS, F_batch)
    model_batch.set_batch_map(F_batch.batch)
    batched_mg, batched_g = CMGDB.ComputeMorseGraph(model_batch)
    assert morse_signature(reference_mg, reference_g) == morse_signature(batched_mg, batched_g)
