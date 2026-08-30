import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

import CMGDB
from CMGDB.PlotMorseSets import _exposed_faces, _grid_index, _merged_outline, _outline_loops


def signed_area(loop):
    x, y = loop[:, 0], loop[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def test_outline_loops_hole_and_corner_touch():
    # A 5x5 ring round a 3x3 hole, plus two cells touching it and each other
    # at a corner only.
    cells = [(i, j) for i in range(5) for j in range(5) if not (1 <= i <= 3 and 1 <= j <= 3)]
    cells += [(5, 5), (6, 6)]
    loops = _outline_loops(np.array(cells))
    assert len(loops) == 4
    # The hole winds clockwise (negative area), everything else counter-clockwise.
    assert sorted(signed_area(l) for l in loops) == pytest.approx([-9.0, 1.0, 1.0, 25.0])
    # Collinear vertices are dropped: every loop here is a plain rectangle.
    assert all(len(l) == 4 for l in loops)


def test_outline_loops_straight_run_is_one_segment():
    loops = _outline_loops(np.array([(k, 0) for k in range(50)]))
    assert len(loops) == 1 and len(loops[0]) == 4


def test_outline_loops_ignore_duplicate_cells():
    loops = _outline_loops(np.array([(0, 0), (0, 0), (1, 0)]))
    assert len(loops) == 1 and signed_area(loops[0]) == pytest.approx(2.0)


def test_grid_index_rejects_mixed_widths():
    lower = np.array([[0.0, 0.0], [1.0, 0.0]])
    upper = np.array([[1.0, 1.0], [1.5, 1.0]])
    assert _grid_index(lower, upper) is None
    assert _merged_outline(lower, upper) is None


H = 0.05


def blob_rows(seed=0):
    """Two Morse sets of aligned boxes: an annulus round a hole, and scattered cells."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(-14, 15):
        for j in range(-14, 15):
            if 5 < np.hypot(i + 0.5, j + 0.5) < 12:
                rows.append([i * H, j * H, (i + 1) * H, (j + 1) * H, 0])
    for _ in range(40):
        i, j = rng.integers(15, 25, size=2)
        rows.append([i * H, j * H, (i + 1) * H, (j + 1) * H, 1])
    return rows


def render(fig):
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].astype(int)


def pixel(image, ax, x, y):
    px, py = ax.transData.transform((x, y))
    return image[int(image.shape[0] - py), int(px)]


def test_merged_plot_matches_per_box_plot():
    rows = blob_rows()
    fig_m, ax_m = CMGDB.PlotMorseSets(rows, show=False)
    fig_b, ax_b = CMGDB.PlotMorseSets(rows, merge_boxes=False, show=False)
    # One path per Morse set against one collection of boxes per Morse set.
    assert len(ax_m.patches) == 2 and len(ax_m.collections) == 0
    assert len(ax_b.collections) == 2 and len(ax_b.patches) == 0
    a, b = render(fig_m), render(fig_b)
    # Only antialiasing along the staircase boundary may differ.
    assert (np.abs(a - b).max(axis=2) > 8).mean() < 0.005
    # The hole is left empty under the nonzero fill rule, and the ring is filled.
    assert (pixel(a, ax_m, 0.0, 0.0) == 255).all()
    assert (pixel(a, ax_m, 8.5 * H, 0.5 * H) != 255).any()
    plt.close(fig_m)
    plt.close(fig_b)


def test_zoom_inset_uses_merged_paths():
    rows = blob_rows()
    fig, ax = CMGDB.PlotMorseSets(rows, zoom_nodes=[1], show=False)
    inset = ax.child_axes[0]
    assert len(inset.patches) == 2 and len(inset.collections) == 0
    plt.close(fig)


@pytest.mark.parametrize("kwargs, merged, per_box", [
    ({'merge_boxes': False}, 0, 2),
    ({'edge_clr': 'k'}, 0, 2),             # outlined boxes are drawn one by one
    ({'scale_factor': [1.5, 1]}, 1, 1),    # an inflated set leaves the grid
])
def test_per_box_fallbacks(kwargs, merged, per_box):
    fig, ax = CMGDB.PlotMorseSets(blob_rows(), show=False, **kwargs)
    assert len(ax.patches) == merged and len(ax.collections) == per_box
    plt.close(fig)


def test_unaligned_set_falls_back_alone():
    rows = blob_rows() + [[2.0, 2.0, 2.0 + 1.5 * H, 2.0 + H, 1]]
    fig, ax = CMGDB.PlotMorseSets(rows, show=False)
    assert len(ax.patches) == 1 and len(ax.collections) == 1
    plt.close(fig)


def test_dpi_defaults_follow_rasterize(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(matplotlib.figure.Figure, 'savefig',
                        lambda self, fname, **kw: seen.append(kw.get('dpi')))
    rows = blob_rows()
    CMGDB.PlotMorseSets(rows, fig_fname=str(tmp_path / 'a.pdf'), show=False)
    assert seen[-1] == 300
    CMGDB.PlotMorseSets(rows, fig_fname=str(tmp_path / 'b.pdf'), rasterize=True, show=False)
    assert seen[-1] == 600
    CMGDB.PlotMorseSets(rows, fig_fname=str(tmp_path / 'c.pdf'), rasterize=True, dpi=200,
                        show=False)
    assert seen[-1] == 200
    plt.close('all')


def cube_rows(n=3):
    return [[i * H, j * H, k * H, (i + 1) * H, (j + 1) * H, (k + 1) * H, 0]
            for i in range(n) for j in range(n) for k in range(n)]


def test_exposed_faces_still_cull_interior():
    faces, labels, aligned = _exposed_faces(cube_rows(), [0], [1])
    assert aligned and len(faces) == 6 * 9


def test_3d_rasterize_is_decided_by_face_count(monkeypatch, tmp_path):
    import importlib
    plot_module = importlib.import_module('CMGDB.PlotMorseSets')
    rows = cube_rows()                       # 54 exposed faces: stays vector
    out = tmp_path / 'vector.pdf'
    fig, ax = CMGDB.PlotMorseSets3D(rows, fig_fname=str(out), show=False)
    assert not ax.collections[0].get_rasterized()
    assert b'/Subtype /Image' not in out.read_bytes()
    plt.close(fig)
    # Above the threshold the faces become one bitmap inside the vector page.
    monkeypatch.setattr(plot_module, 'RASTERIZE_FACES', 40)
    out = tmp_path / 'raster.pdf'
    fig, ax = CMGDB.PlotMorseSets3D(rows, fig_fname=str(out), show=False)
    assert ax.collections[0].get_rasterized()
    assert b'/Subtype /Image' in out.read_bytes()
    plt.close(fig)
    # An explicit value overrides the count either way.
    fig, ax = CMGDB.PlotMorseSets3D(rows, rasterize=False, show=False)
    assert not ax.collections[0].get_rasterized()
    plt.close(fig)
    monkeypatch.setattr(plot_module, 'RASTERIZE_FACES', 100000)
    fig, ax = CMGDB.PlotMorseSets3D(rows, rasterize=True, show=False)
    assert ax.collections[0].get_rasterized()
    plt.close(fig)
