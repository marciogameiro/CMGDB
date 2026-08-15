"""Tests for PrecomputedBoxMap: lookups must reproduce live BoxMap evaluation
exactly (corner mode) on every dyadic box of the subdivision tree, and the
full Morse graph computed through the precomputed table must match the one
computed with live batched evaluation."""

import numpy as np
import pytest
import CMGDB

THETA = np.array([19.6, 23.68])
LOWER_BOUNDS = [-0.001, -0.001]
UPPER_BOUNDS = [90.0, 70.0]


def f_vec(X):
    s = X.sum(axis=1)
    y0 = (X @ THETA) * np.exp(-0.1 * s)
    return np.column_stack([y0, 0.7 * X[:, 0]])


def f_scalar(x):
    return list(f_vec(np.array([x]))[0])


# Map using only IEEE-exact operations (+, -, *, /): unlike np.exp -- whose
# SIMD (large batch) and scalar (single point) code paths may differ in the
# last ulp on some platforms -- this map is bit-identical for any evaluation
# batching, so lookups can be asserted strictly equal to live evaluation.
def g_vec(X):
    return X * (90.0 - X) / 45.0


def g_scalar(x):
    return list(g_vec(np.array([x]))[0])


SUBDIV_MAX = 12


@pytest.fixture(scope="module")
def F_pre():
    return CMGDB.PrecomputedBoxMap(f_vec, LOWER_BOUNDS, UPPER_BOUNDS, SUBDIV_MAX)


@pytest.fixture(scope="module")
def G_pre():
    return CMGDB.PrecomputedBoxMap(g_vec, LOWER_BOUNDS, UPPER_BOUNDS, SUBDIV_MAX)


def dyadic_boxes(F, rng, count):
    """Random boxes aligned to the subdiv_max lattice, at assorted depths,
    constructed with the same arithmetic as the lattice coordinates."""
    cells = F._cells_per_axis
    side = F._finest_box_side
    lower = np.asarray(LOWER_BOUNDS, dtype=float)
    boxes = []
    for _ in range(count):
        spans = np.array([cells[d] >> rng.integers(0, int(cells[d]).bit_length())
                          for d in range(2)], dtype=np.int64)
        spans = np.maximum(spans, 1)
        i_lower = np.array([rng.integers(0, cells[d] // spans[d]) * spans[d]
                            for d in range(2)], dtype=np.int64)
        i_upper = i_lower + spans
        rect = list(lower + i_lower * side) + list(lower + i_upper * side)
        boxes.append(rect)
    return boxes


def test_corner_mode_exact_vs_live_boxmap(G_pre):
    rng = np.random.default_rng(23)
    for rect in dyadic_boxes(G_pre, rng, 200):
        assert G_pre(rect) == CMGDB.BoxMap(g_scalar, rect)


def test_corner_mode_with_padding_exact(G_pre):
    G_pad = CMGDB.PrecomputedBoxMap(g_vec, LOWER_BOUNDS, UPPER_BOUNDS, SUBDIV_MAX,
                                    padding=True)
    rng = np.random.default_rng(29)
    for rect in dyadic_boxes(G_pad, rng, 100):
        assert G_pad(rect) == CMGDB.BoxMap(g_scalar, rect, padding=True)


def test_corner_mode_close_for_exp_map(F_pre):
    # np.exp evaluated in a large batch may differ from single-point
    # evaluation in the last ulp (SIMD vs scalar code paths), so for maps
    # using transcendental functions the guarantee is agreement up to
    # floating-point evaluation-batching effects
    rng = np.random.default_rng(23)
    for rect in dyadic_boxes(F_pre, rng, 200):
        assert np.allclose(F_pre(rect), CMGDB.BoxMap(f_scalar, rect),
                           rtol=1e-12, atol=1e-12)


def test_center_mode_matches_live(F_pre):
    F_center = CMGDB.PrecomputedBoxMap(f_vec, LOWER_BOUNDS, UPPER_BOUNDS, SUBDIV_MAX,
                                       mode='center')
    rng = np.random.default_rng(31)
    for rect in dyadic_boxes(F_center, rng, 100):
        expected = CMGDB.BoxMap(f_scalar, rect, mode='center')
        assert np.allclose(F_center(rect), expected, rtol=1e-12, atol=1e-12)


def test_whole_domain_box(G_pre):
    rect = list(LOWER_BOUNDS) + list(UPPER_BOUNDS)
    assert G_pre(rect) == CMGDB.BoxMap(g_scalar, rect)


def test_batch_matches_scalar(F_pre):
    rng = np.random.default_rng(37)
    rects = dyadic_boxes(F_pre, rng, 50)
    batched = F_pre.batch(rects)
    for i, rect in enumerate(rects):
        assert list(batched[i]) == F_pre(rect)


def test_chunking_identical():
    # Chunk boundaries must not change the table (exact-arithmetic map, so
    # SIMD batching effects cannot mask a real chunking bug)
    small_chunks = CMGDB.PrecomputedBoxMap(g_vec, LOWER_BOUNDS, UPPER_BOUNDS, 8,
                                           batch_points=100)
    auto_chunks = CMGDB.PrecomputedBoxMap(g_vec, LOWER_BOUNDS, UPPER_BOUNDS, 8)
    assert np.array_equal(small_chunks._table, auto_chunks._table)


def test_box_finer_than_lattice_raises(F_pre):
    side = F_pre._finest_box_side
    rect = [0.0, 0.0, side[0] / 2, side[1] / 2]
    with pytest.raises(ValueError):
        F_pre(rect)


def test_off_lattice_box_raises(F_pre):
    side = F_pre._finest_box_side
    rect = [0.3 * side[0], 0.0, 2.3 * side[0], side[1]]
    with pytest.raises(ValueError):
        F_pre(rect)


def test_invalid_mode_raises():
    with pytest.raises(ValueError):
        CMGDB.PrecomputedBoxMap(f_vec, LOWER_BOUNDS, UPPER_BOUNDS, 8, mode='random')


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


def test_morse_graph_matches_live_evaluation():
    subdiv_min, subdiv_max = 14, 16

    def F_live(rect):
        return CMGDB.BoxMap(f_scalar, rect)

    def F_live_batch(rects):
        return CMGDB.BoxMapBatch(f_vec, rects)

    model = CMGDB.Model(subdiv_min, subdiv_max, LOWER_BOUNDS, UPPER_BOUNDS, F_live)
    model.set_batch_map(F_live_batch)
    live_mg, live_g = CMGDB.ComputeMorseGraph(model)

    F = CMGDB.PrecomputedBoxMap(f_vec, LOWER_BOUNDS, UPPER_BOUNDS, subdiv_max)
    model_pre = CMGDB.Model(subdiv_min, subdiv_max, LOWER_BOUNDS, UPPER_BOUNDS, F)
    model_pre.set_batch_map(F.batch)
    pre_mg, pre_g = CMGDB.ComputeMorseGraph(model_pre)

    assert morse_signature(live_mg, live_g) == morse_signature(pre_mg, pre_g)


def test_torch_module_path():
    torch = pytest.importorskip("torch")

    class Net(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.linear = torch.nn.Linear(2, 2)

        def forward(self, x):
            return torch.sigmoid(self.linear(x)) * 50.0

    torch.manual_seed(0)
    net = Net()

    def f_numpy(P):
        with torch.no_grad():
            x = torch.as_tensor(P, dtype=torch.float32)
            return (torch.sigmoid(net.linear(x)) * 50.0).numpy().astype(float)

    F_torch = CMGDB.PrecomputedBoxMap(net, LOWER_BOUNDS, UPPER_BOUNDS, 8, device='cpu')
    F_ref = CMGDB.PrecomputedBoxMap(f_numpy, LOWER_BOUNDS, UPPER_BOUNDS, 8)
    assert np.array_equal(F_torch._table, F_ref._table)
