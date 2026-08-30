### PlotMorseSets.py
### MIT LICENSE 2026 Marcio Gameiro

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path
import CMGDB

# Default color list
DEFAULT_CLIST = ['#1f77b4', '#e6550d', '#31a354', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                 '#bcbd22', '#80b1d3', '#ffffb3', '#fccde5', '#b3de69', '#fdae6b', '#6a3d9a', '#c49c94',
                 '#fb8072', '#dbdb8d', '#bc80bd', '#ffed6f', '#637939', '#c5b0d5', '#636363', '#c7c7c7',
                 '#8dd3c7', '#b15928', '#e8cb32', '#9e9ac8', '#74c476', '#ff7f0e', '#9edae5', '#90d743',
                 '#e7969c', '#17becf', '#7b4173', '#8ca252', '#ad494a', '#8c6d31', '#a55194', '#00cc49']

# Face count above which PlotMorseSets3D rasterizes by default. Measured on
# cubical tubes: vector output costs 24-30 bytes a face, the 600 dpi bitmap
# 1.4 MB at 25k faces and 2.1 MB at 75k, where the two are the same size. Vector
# is preferred at a tie, having no resolution ceiling, so the switch sits above.
RASTERIZE_FACES = 100000


def _load_morse_sets(morse_sets):
    """Accept a Morse graph, a saved Morse set file name, or a list of boxes.

    Returns the box rows and the number of Morse sets, or None for the count
    when it can only be inferred from the labels present.
    """
    if type(morse_sets) == CMGDB._cmgdb.MorseGraph:
        morse_graph = morse_sets
        num_morse_sets = morse_graph.num_vertices()
        rows = [box + [node] for node in range(num_morse_sets)
                for box in morse_graph.morse_set_boxes(node)]
        return rows, num_morse_sets
    if type(morse_sets) == str:
        return CMGDB.LoadMorseSetFile(morse_sets), None
    return list(morse_sets), None


def _box_dim(rows):
    rect = rows[0]
    assert len(rect) % 2 == 1, "Wrong dimension in Morse sets data"
    return int((len(rect) - 1) / 2)


def _resolve_plot_setup(rows, num_morse_sets, morse_nodes, cmap, clist, scale_factor):
    """Shared colour map, node list, and per-node scale factors."""
    if num_morse_sets == None:
        num_morse_sets = max([int(rect[-1]) for rect in rows]) + 1
    if morse_nodes == None:
        morse_nodes = range(num_morse_sets)
    if scale_factor == None:
        scale_factor = [1] * num_morse_sets
    elif len(scale_factor) < num_morse_sets:
        # Indexing is by Morse node, so a short list would raise on the first
        # unlisted set. Pad rather than fail: unlisted sets are drawn faithfully.
        scale_factor = list(scale_factor) + [1] * (num_morse_sets - len(scale_factor))
    if cmap == None and clist == None:
        clist = DEFAULT_CLIST
    if cmap == None:
        cmap = matplotlib.colors.ListedColormap(clist[:num_morse_sets])
    try:
        num_colors = len(cmap.colors)
    except:
        num_colors = 0
    if (num_colors > 0) and (num_colors < num_morse_sets):
        cmap_norm = lambda k: k % num_colors
    else:
        cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=num_morse_sets - 1)
    return num_morse_sets, morse_nodes, scale_factor, cmap, cmap_norm



def _effective_dpi_scale(fig, dpi, probe_dpi=50):
    """Factor that turns a requested dpi into the dpi seen on the saved page.

       Matplotlib rasterizes an artist over the region it draws, while
       ``bbox_inches='tight'`` crops the page to the content plus its labels.
       The bitmap covers fewer inches than the page, so its resolution there is
       lower than the number requested -- by about 17% for a 2-D plot and 49%
       for a 3-D one, where the projected axes fill a fraction of the canvas.

       The ratio depends on the layout, not on the dpi, so one cheap probe
       render measures it exactly. Modelling it from the axes or the projected
       data cube does not work: neither predicts what mplot3d actually emits.
    """
    import io
    import re

    try:
        buffer = io.BytesIO()
        fig.savefig(buffer, format='pdf', dpi=probe_dpi, bbox_inches='tight')
        raw = buffer.getvalue()
        widths = [int(m) for m in re.findall(rb'/Width\s+(\d+)', raw)]
        boxes = re.findall(rb'/MediaBox\s*\[([^\]]+)\]', raw)
        if not widths or not boxes:
            return 1.0
        numbers = [float(v) for v in boxes[0].split()]
        page_in = (numbers[2] - numbers[0]) / 72.0
        bitmap_in = max(widths) / probe_dpi
        if bitmap_in <= 0 or page_in <= 0:
            return 1.0
        return max(1.0, page_in / bitmap_in)
    except Exception:
        # A backend that cannot probe gets the uncorrected dpi rather than an
        # exception raised from a plotting call.
        return 1.0


def _finish(fig, ax, fig_fname, dpi, show, rasterized=False):
    """Save if asked, show if the backend can, and hand the figure back.

    ``show=None`` shows only on an interactive backend. Under Agg -- scripts,
    CI, batch rendering -- an unconditional ``plt.show()`` does nothing except
    emit a warning on every call, while in a notebook the figure appears
    anyway. Returning ``(fig, ax)`` lets a caller adjust the plot before
    saving, which was previously impossible.

    ``dpi=None`` resolves to 600 when the plot is rasterized and to 300
    otherwise. On a vector page the number only sets the resolution of the
    embedded bitmap, where 600 keeps it sharp in print at little cost; a
    fully vector page ignores it, and a bitmap format keeps its former size.
    """
    if dpi == None:
        dpi = 600 if rasterized else 300
    if fig_fname:
        # Scale so ``dpi`` is the resolution obtained on the page, not the one
        # requested of a bitmap that covers only part of it.
        save_dpi = dpi * (_effective_dpi_scale(fig, dpi) if rasterized else 1.0)
        fig.savefig(fig_fname, dpi=save_dpi, bbox_inches='tight')
    if show == None:
        show = matplotlib.get_backend().lower() not in ('agg', 'pdf', 'ps', 'svg', 'cairo', 'template')
    if show:
        plt.show()
    return fig, ax


def _drawn_span(rows, morse_nodes, dim, d, scale_factor):
    """Extent actually drawn along axis d, honouring per-set inflation.

       scale_factor enlarges a box about its own centre, so a scaled box
       reaches beyond its true endpoints. Limits taken from the raw endpoints
       would clip exactly the small sets that were inflated to be visible.
    """
    lo, hi = None, None
    for rect in rows:
        node = int(rect[-1])
        if node not in morse_nodes:
            continue
        a, b = rect[d], rect[dim + d]
        if scale_factor != None:
            factor = scale_factor[node] if node in scale_factor else 1.0
            centre, half = (a + b) / 2, (b - a) * factor / 2
            a, b = centre - half, centre + half
        lo = a if lo == None or a < lo else lo
        hi = b if hi == None or b > hi else hi
    return lo, hi


def _axis_limits(rows, morse_nodes, dim, d1, d2, xlim, ylim,
                 scale_factor=None, margin=0.0):
    """Axis limits for a two-dimensional Morse set plot.

       margin is a fraction of the drawn span added to each side, so the sets
       sit inside the axes instead of running into the bounding box. Explicit
       xlim/ylim are used exactly as given -- a caller that names its limits
       means them.
    """
    if xlim == None:
        x_min, x_max = _drawn_span(rows, morse_nodes, dim, d1, scale_factor)
        if x_max - x_min < 0.1:
            x_min -= 0.05
            x_max += 0.05
        pad = (x_max - x_min) * margin
        x_min, x_max = x_min - pad, x_max + pad
    else:
        x_min, x_max = xlim[0], xlim[1]
    if ylim == None:
        y_min, y_max = _drawn_span(rows, morse_nodes, dim, d2, scale_factor)
        if y_max - y_min < 0.1:
            y_min -= 0.05
            y_max += 0.05
        pad = (y_max - y_min) * margin
        y_min, y_max = y_min - pad, y_max + pad
    else:
        y_min, y_max = ylim[0], ylim[1]
    return x_min, x_max, y_min, y_max


def _grid_index(lower, upper):
    """Integer grid coordinates of boxes that share one aligned grid.

       CMGDB's subdivision produces boxes of one size on one lattice, which is
       what lets adjacency be decided by index arithmetic instead of geometry.
       Returns (origin, cell, index) with lower == origin + index * cell, or
       None when the boxes do not fit one grid -- mixed depths, say -- and the
       caller must fall back to per-box geometry.
    """
    widths = upper - lower
    cell = np.median(widths, axis=0)
    if not (np.all(cell > 0) and np.allclose(widths, cell, rtol=1e-8, atol=1e-11)):
        return None
    origin = lower.min(axis=0)
    index = np.rint((lower - origin) / cell).astype(np.int64)
    if not np.allclose(lower, origin + index * cell, rtol=1e-8, atol=1e-10):
        return None
    return origin, cell, index


def _outline_loops(index):
    """Boundary of a set of grid cells as closed loops of vertex indices.

       Cell (i, j) covers [i, i+1] x [j, j+1]. A side lies on the boundary
       exactly when the cell across it is absent, so the boundary is found by
       index lookups: the two-dimensional counterpart of _exposed_faces. Each
       boundary side is directed with its own cell on the left; chained head
       to tail, outer boundaries then wind counter-clockwise and holes
       clockwise, which is what the nonzero fill rule of the PDF and Agg
       backends needs to leave a hole empty. Where two components touch at a
       corner the left turn keeps each on its own loop; every side is used
       exactly once with that winding, so a loop that does pass through such a
       vertex still fills correctly. Collinear vertices are dropped, so a
       straight run of cells costs one segment however long it is.
    """
    index = np.unique(np.asarray(index, dtype=np.int64), axis=0)
    if len(index) == 0:
        return []
    i0, j0 = index.min(axis=0)
    i, j = index[:, 0] - i0, index[:, 1] - j0
    stride = j.max() + 2            # j + 1 stays below it, so keys never collide
    occupied = np.sort(i * stride + j)
    sides = []
    # (neighbour offset, side start, side end): each side runs so that its
    # own cell lies on the left
    for (di, dj), (ax, ay), (bx, by) in (((0, -1), (0, 0), (1, 0)),    # bottom, +x
                                         ((1, 0), (1, 0), (1, 1)),     # right, +y
                                         ((0, 1), (1, 1), (0, 1)),     # top, -x
                                         ((-1, 0), (0, 1), (0, 0))):   # left, -y
        exposed = ~np.isin((i + di) * stride + (j + dj), occupied)
        sides.append(np.column_stack([i[exposed] + ax, j[exposed] + ay,
                                      i[exposed] + bx, j[exposed] + by]))
    sides = np.concatenate(sides)

    vstride = j.max() + 3           # vertex coordinates run one past the cells
    leaving = {}
    for e, k in enumerate((sides[:, 0] * vstride + sides[:, 1]).tolist()):
        leaving.setdefault(k, []).append(e)
    used = np.zeros(len(sides), dtype=bool)
    loops = []
    for first in range(len(sides)):
        if used[first]:
            continue
        loop, e = [], first
        while True:
            used[e] = True
            loop.append(sides[e, :2])
            end = sides[e, 2:]
            candidates = [c for c in leaving.get(int(end[0] * vstride + end[1]), [])
                          if not used[c]]
            if not candidates:
                break               # back at the start: in and out degrees match
            if len(candidates) == 1:
                e = candidates[0]
            else:
                # Two components meet at this corner. The interior is on the
                # left, so the left turn continues along this component.
                step = end - sides[e, :2]
                left = np.array([-step[1], step[0]])
                e = max(candidates, key=lambda c: int((sides[c, 2:] - sides[c, :2]) @ left))
        loop = np.array(loop)
        # A vertex between two sides of the same direction is not a corner.
        step = np.roll(loop, -1, axis=0) - loop
        corner = np.any(step != np.roll(step, 1, axis=0), axis=1)
        loops.append(loop[corner] + (i0, j0))
    return loops


def _merged_outline(lower, upper):
    """One compound path covering the union of aligned boxes, or None.

       Drawing a Morse set as the outline of its boxes rather than as the boxes
       themselves is what keeps the file small: a set of 10^5 cells is a
       handful of polygons, and every interior edge was invisible anyway.
    """
    grid = _grid_index(lower, upper)
    if grid is None:
        return None
    origin, cell, index = grid
    vertices, codes = [], []
    for loop in _outline_loops(index):
        points = origin + loop * cell
        vertices.extend(points.tolist())
        vertices.append(points[0].tolist())
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(points) - 1) + [Path.CLOSEPOLY])
    return Path(vertices, codes)


def _draw_boxes(ax, rows, morse_nodes, dim, d1, d2, cmap, cmap_norm,
                scale_factor, edge_clr, linewidth, alpha, rasterize, merge_boxes=True):
    """Add one artist per Morse set to ax. Returns the artists.

       By default a set is one PathPatch, the outline of the union of its
       boxes (see _merged_outline). The boxes are drawn one by one, as a
       PatchCollection, when merging is switched off, when the boxes do not
       share one grid, when a scale_factor other than 1 moves them off it, or
       when an explicit edge_clr asks for each box to be outlined.
    """
    drawn = []
    for morse_node in morse_nodes:
        morse_set = [rect for rect in rows if int(rect[-1]) == morse_node]
        if not morse_set:
            continue
        clr = matplotlib.colors.to_hex(cmap(cmap_norm(morse_node)), keep_alpha=True)
        factor = scale_factor[morse_node]
        # Edge in the face colour by default: it closes the antialiasing seam
        # between neighbouring boxes instead of outlining each one.
        edges = clr if edge_clr == None else edge_clr
        path = None
        if merge_boxes and edge_clr == None and factor == 1:
            boxes = np.asarray([[float(v) for v in rect] for rect in morse_set])
            path = _merged_outline(boxes[:, [d1, d2]], boxes[:, [dim + d1, dim + d2]])
        if path is not None:
            artist = PathPatch(path, facecolor=clr, edgecolor=edges, linewidth=linewidth,
                               alpha=alpha, rasterized=rasterize)
            ax.add_patch(artist)
        else:
            patches = []
            for rect in morse_set:
                x0, y0 = rect[d1], rect[d2]
                x1, y1 = rect[dim + d1], rect[dim + d2]
                width, height = (x1 - x0) * factor, (y1 - y0) * factor
                cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
                patches.append(Rectangle((cx - width / 2, cy - height / 2), width, height))
            # One collection per Morse set rather than one artist per box: the
            # box count reaches 10^5 routinely and per-artist overhead
            # dominates there.
            artist = PatchCollection(patches, facecolors=clr, edgecolors=edges,
                                     linewidths=linewidth, alpha=alpha,
                                     rasterized=rasterize)
            ax.add_collection(artist)
        drawn.append(artist)
    return drawn


def _zoom_region(rows, zoom_nodes, dim, d1, d2, scale_factor, pad, square):
    """Data-space window enclosing zoom_nodes, padded and optionally squared."""
    xs = [r for r in rows if int(r[-1]) in zoom_nodes]
    if not xs:
        raise ValueError("zoom_nodes match no boxes: " + repr(sorted(zoom_nodes)))
    x0 = min(r[d1] for r in xs)
    x1 = max(r[dim + d1] for r in xs)
    y0 = min(r[d2] for r in xs)
    y1 = max(r[dim + d2] for r in xs)
    if square:
        # A Morse set can be a near-degenerate sliver; a square window keeps the
        # inset from stretching one axis by orders of magnitude against the other.
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        half = max(x1 - x0, y1 - y0) / 2
        x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half
    px, py = (x1 - x0) * pad, (y1 - y0) * pad
    return x0 - px, x1 + px, y0 - py, y1 + py


def _inset_corner(ax, region):
    """Corner of the axes furthest from the zoomed window, in axes fraction."""
    (xa, xb), (ya, yb) = ax.get_xlim(), ax.get_ylim()
    fx = ((region[0] + region[1]) / 2 - xa) / (xb - xa) if xb != xa else 0.5
    fy = ((region[2] + region[3]) / 2 - ya) / (yb - ya) if yb != ya else 0.5
    size = 0.38
    left = 0.04 if fx > 0.5 else 1.0 - size - 0.04
    bottom = 0.04 if fy > 0.5 else 1.0 - size - 0.04
    return [left, bottom, size, size]


def _add_zoom_inset(ax, rows, morse_nodes, dim, d1, d2, cmap, cmap_norm,
                    scale_factor, edge_clr, linewidth, alpha, rasterize, merge_boxes,
                    zoom_nodes, zoom_bounds, zoom_pos, zoom_pad, zoom_square,
                    zoom_edge_clr, zoom_lw, zoom_ticks, fontsize):
    """Draw a magnified copy of one region, boxed and joined to its source.

       Morse sets differ in extent by orders of magnitude, so a set that is
       only a few boxes across is invisible next to a large one. Inflating it
       with scale_factor makes it visible but draws it at the wrong size; an
       inset instead shows it at true relative geometry, magnified, with the
       box and connector lines saying where it came from.
    """
    if zoom_bounds != None:
        region = tuple(zoom_bounds)
    else:
        nodes = set(int(n) for n in zoom_nodes)
        region = _zoom_region(rows, nodes, dim, d1, d2, scale_factor,
                              zoom_pad, zoom_square)
    position = _inset_corner(ax, region) if zoom_pos == None else list(zoom_pos)
    axins = ax.inset_axes(position)
    # Every set is redrawn, not just the zoomed ones: neighbouring structure is
    # what makes the magnified view legible as part of the whole picture.
    _draw_boxes(axins, rows, morse_nodes, dim, d1, d2, cmap, cmap_norm,
                scale_factor, edge_clr, linewidth, alpha, rasterize, merge_boxes)
    axins.set_xlim(region[0], region[1])
    axins.set_ylim(region[2], region[3])
    if zoom_ticks:
        axins.tick_params(labelsize=max(fontsize - 5, 5))
    else:
        axins.set_xticks([])
        axins.set_yticks([])
    for spine in axins.spines.values():
        spine.set_edgecolor(zoom_edge_clr)
        spine.set_linewidth(zoom_lw)
    ax.indicate_inset_zoom(axins, edgecolor=zoom_edge_clr, linewidth=zoom_lw,
                           alpha=1.0)
    return axins


def PlotMorseSets(morse_sets, morse_nodes=None, proj_dims=None, cmap=None, clist=None,
                  scale_factor=None, fig_w=8, fig_h=8, xlim=None, ylim=None, margin=0.02,
                  axis_labels=True, xlabel='$x$', ylabel='$y$', fontsize=15, edge_clr=None,
                  linewidth=0.5, alpha=None, zoom_nodes=None, zoom_bounds=None,
                  zoom_pos=None, zoom_pad=0.25, zoom_square=True, zoom_edge_clr='0.35',
                  zoom_lw=0.8, zoom_ticks=False, merge_boxes=True, fig_fname=None,
                  dpi=None, rasterize=False, show=None):
    """Plot Morse sets as filled rectangles, one per box.

       Each box is drawn at its true extent, so the picture is exact at any
       zoom: unlike a scatter of square markers, a rectangle keeps a box's
       aspect ratio and does not overlap its neighbours. This is the default
       Morse set plot; PlotMorseSetsScatter draws the older marker version,
       which stays useful when boxes are so small that markers read better.

       margin pads the axes by a fraction of the drawn span on each side so
       the sets do not run into the bounding box; limits follow the scaled
       box extents too, so an inflated set is never clipped. Pass xlim/ylim
       to override the computed limits entirely.

       zoom_nodes magnifies the region holding those Morse sets in an inset,
       with a box drawn round the region and lines joining it to the inset.
       Prefer it to a large scale_factor when a set is too small to see: the
       inset keeps every box at its true size and shape, whereas inflation
       draws the set larger than it is. zoom_bounds names the window directly
       as (x0, x1, y0, y1); zoom_pos places the inset as [x0, y0, w, h] in axes
       fractions, defaulting to the corner furthest from the region.

       morse_sets is a Morse graph, a saved Morse set file name, or a list of
       boxes [lower..., upper..., label]. scale_factor is a list indexed by
       Morse node that enlarges a set about each box's centre, so a set orders
       of magnitude smaller than another stays visible; positions and relative
       geometry are unchanged, only the drawn size.

       edge_clr colours the rectangle borders. The default draws them in each
       box's own face colour, so adjacent boxes of one Morse set read as a
       single region: with no edge at all, antialiasing leaves pale seams
       between neighbours. Give it a colour to outline the boxes instead.

       alpha applies to the whole collection, faces and edges together, so a
       translucent Morse set does not darken along the seams where an edge is
       drawn over its own face. Leave it None to use the alpha carried by the
       colours themselves: an 8-digit hex in clist such as '#1f77b480' is
       honoured as given, and needs no unpacking. A scalar alpha overrides any
       per-colour alpha, which is why it is not set by default.

       merge_boxes draws each Morse set as a single path, the outline of the
       union of its boxes, rather than as one rectangle per box. The picture
       is the same -- the interior edges were never visible -- but a set of
       10^5 boxes becomes a few polygons, so the file is orders of magnitude
       smaller and opens at once. Set it False to draw the boxes one by one
       as before; that also happens by itself when the boxes do not share one
       grid, when a scale_factor other than 1 moves them off it, or when an
       explicit edge_clr asks for each box to be outlined.

       rasterize draws the boxes as a raster image inside the vector figure,
       at dpi (600 by default when rasterizing). With merge_boxes the vector
       output is already small, so this mainly serves the fallback cases
       above; it costs a resolution ceiling.

       Returns (fig, ax).
    """
    rows, num_morse_sets = _load_morse_sets(morse_sets)
    dim = _box_dim(rows)
    if dim == 1:
        return PlotMorseSets1D(morse_sets, morse_nodes=morse_nodes, cmap=cmap, clist=clist,
                               fig_w=fig_w, xlim=xlim, axis_labels=axis_labels,
                               xlabel=xlabel, fontsize=fontsize, alpha=alpha,
                               fig_fname=fig_fname, dpi=dpi, rasterize=rasterize, show=show)
    num_morse_sets, morse_nodes, scale_factor, cmap, cmap_norm = _resolve_plot_setup(
        rows, num_morse_sets, morse_nodes, cmap, clist, scale_factor)
    if proj_dims == None:
        d1, d2 = 0, 1
    else:
        d1, d2 = proj_dims[0], proj_dims[1]
    assert max(d1, d2) < dim, "Wrong projection dimensions"

    x_min, x_max, y_min, y_max = _axis_limits(rows, morse_nodes, dim, d1, d2, xlim, ylim,
                                              scale_factor=scale_factor, margin=margin)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    drawn = _draw_boxes(ax, rows, morse_nodes, dim, d1, d2, cmap, cmap_norm,
                        scale_factor, edge_clr, linewidth, alpha, rasterize, merge_boxes)

    if zoom_nodes != None or zoom_bounds != None:
        _add_zoom_inset(ax, rows, morse_nodes, dim, d1, d2, cmap, cmap_norm,
                        scale_factor, edge_clr, linewidth, alpha, rasterize, merge_boxes,
                        zoom_nodes, zoom_bounds, zoom_pos, zoom_pad, zoom_square,
                        zoom_edge_clr, zoom_lw, zoom_ticks, fontsize)

    if axis_labels:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    return _finish(fig, ax, fig_fname, dpi, show, rasterize)


def PlotMorseSets1D(morse_sets, morse_nodes=None, cmap=None, clist=None, scale_factor=None,
                    fig_w=8, fig_h=None, xlim=None, axis_labels=True, xlabel='$x$',
                    axis_arrow=True, fontsize=15, height=0.18, edge_clr=None,
                    linewidth=0.5, alpha=None, label_sets=True, merge_tol=1e-9,
                    fig_fname=None, dpi=None, rasterize=False, show=None):
    """Plot one-dimensional Morse sets as boxes sitting on the axis.

       Every set is drawn on the same line, which is what a one-dimensional
       decomposition is: disjoint intervals of one space. Boxes straddle the
       axis and are drawn behind it, so the axis line runs through them and
       reads as the phase space itself rather than as a baseline they sit on.

       height is the box height in data units. The axis line sits at the middle
       of the boxes, so the lower half would otherwise run over the tick
       numbers printed beneath it; the tick padding is set from height to keep
       the numbers clear whatever height is chosen. The default is small enough
       that the boxes read as marks on the axis rather than as a filled band.

       Boxes are drawn as rectangles at their true endpoints. A thick line with
       projecting caps would overshoot each interval by half its linewidth in
       display units, which makes short intervals look longer than they are and
       makes the error depend on the figure size.

       axis_arrow draws the axis as an axis: an arrow head at the right end of
       the line, with the label at the tip rather than centred beneath. Tick
       marks are kept either way. Set it False for a plain spine with a
       centred label.

       A Morse set is often several disjoint pieces. Touching boxes are merged
       into maximal intervals -- within merge_tol, since CMGDB's grid gives
       exactly abutting endpoints -- so each piece is one rectangle and each is
       labelled, rather than labelling a midpoint that may fall in a gap.

       Returns (fig, ax).
    """
    rows, num_morse_sets = _load_morse_sets(morse_sets)
    dim = _box_dim(rows)
    assert dim == 1, f"PlotMorseSets1D needs 1-D Morse boxes; got dim={dim}"
    num_morse_sets, morse_nodes, scale_factor, cmap, cmap_norm = _resolve_plot_setup(
        rows, num_morse_sets, morse_nodes, cmap, clist, scale_factor)

    if fig_h == None:
        fig_h = max(1.4, 0.22 * fig_w)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)

    x_lo, x_hi = None, None
    for morse_node in morse_nodes:
        intervals = sorted((rect[0], rect[1]) for rect in rows if int(rect[-1]) == morse_node)
        if not intervals:
            continue
        factor = scale_factor[morse_node]
        # Merge boxes that touch into maximal pieces, so each disjoint piece of
        # the Morse set is one rectangle of the correct length.
        pieces = [list(intervals[0])]
        for a, b in intervals[1:]:
            if a <= pieces[-1][1] + merge_tol:
                pieces[-1][1] = max(pieces[-1][1], b)
            else:
                pieces.append([a, b])
        clr = matplotlib.colors.to_hex(cmap(cmap_norm(morse_node)), keep_alpha=True)
        edges = clr if edge_clr == None else edge_clr
        patches = []
        for a, b in pieces:
            if factor != 1:
                centre, half = (a + b) / 2, (b - a) * factor / 2
                a, b = centre - half, centre + half
            patches.append(Rectangle((a, -height / 2), b - a, height))
            x_lo = a if x_lo == None else min(x_lo, a)
            x_hi = b if x_hi == None else max(x_hi, b)
            if label_sets:
                ax.text((a + b) / 2, height / 2 + 0.06, f'{morse_node}', ha='center',
                        va='bottom', color='black', fontsize=fontsize)
        # zorder below the spines and ticks so the axis line runs over the
        # boxes rather than being hidden by them.
        ax.add_collection(PatchCollection(patches, facecolors=clr, edgecolors=edges,
                                          linewidths=linewidth, alpha=alpha,
                                          rasterized=rasterize, zorder=0))

    if xlim == None:
        pad = 0.02 * (x_hi - x_lo) if x_hi > x_lo else 0.05
        ax.set_xlim([x_lo - pad, x_hi + pad])
    else:
        ax.set_xlim(list(xlim))
    # Symmetric about the axis, with headroom for the labels when drawn.
    ax.set_ylim([-0.5, 0.5 if not label_sets else 0.5 + 0.28])
    ax.set_yticks([])
    # Keep only the bottom spine and put it on the axis line the boxes straddle.
    for spine in ('top', 'left', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['bottom'].set_zorder(3)
    if axis_arrow:
        # Arrow head at the right end of the spine. The marker's transform
        # takes x in axes fractions and y in data units, so it sits on the
        # axis line (y = 0) at the axes edge whatever the data limits are.
        ax.plot(1, 0, '>k', markersize=7, transform=ax.get_yaxis_transform(),
                clip_on=False, zorder=4)
        if axis_labels:
            ax.annotate(xlabel, xy=(1, 0), xycoords=ax.get_yaxis_transform(),
                        xytext=(12, 0), textcoords='offset points', ha='left',
                        va='center', fontsize=fontsize, annotation_clip=False)
    elif axis_labels:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    # The boxes straddle the axis, so their lower half covers the space the
    # tick numbers would occupy. Push the numbers below the boxes, converting
    # the half-height from data units to the points tick_params expects.
    fig.canvas.draw()
    axes_height_points = ax.get_window_extent().height * 72.0 / fig.dpi
    y_lo, y_hi = ax.get_ylim()
    half_box_points = (height / 2.0) / (y_hi - y_lo) * axes_height_points
    ax.tick_params(labelsize=fontsize, pad=half_box_points + 4.0)
    return _finish(fig, ax, fig_fname, dpi, show, rasterize)


def _exposed_faces(rows, morse_nodes, scale_factor):
    """Quadrilateral faces of 3-D boxes that are not shared with a neighbour.

       A solid block of cells is drawn as its shell: interior faces are hidden
       by definition, so emitting them costs time and file size and changes
       nothing on screen. Culling needs the boxes to sit on one aligned grid,
       which is what CMGDB's subdivision produces; when they do not -- mixed
       depths, say -- every face is emitted instead, which is correct but
       larger.

       Faces are scaled about their own box's centre after culling, so a scaled
       set stays a closed surface rather than separating into shells.
    """
    values = np.asarray([[float(v) for v in rect] for rect in rows], dtype=float)
    lower, upper, labels = values[:, :3], values[:, 3:6], values[:, 6].astype(int)
    grid = _grid_index(lower, upper)
    aligned = grid is not None
    occupied = set()
    if aligned:
        index = grid[2]
        occupied = {(int(l), int(i[0]), int(i[1]), int(i[2]))
                    for l, i in zip(labels, index)}
    faces, face_labels = [], []
    for k in range(values.shape[0]):
        label = int(labels[k])
        if label not in morse_nodes:
            continue
        x0, y0, z0 = lower[k]
        x1, y1, z1 = upper[k]
        candidates = (
            ([[x0, y0, z0], [x0, y1, z0], [x0, y1, z1], [x0, y0, z1]], (-1, 0, 0)),
            ([[x1, y0, z0], [x1, y0, z1], [x1, y1, z1], [x1, y1, z0]], (1, 0, 0)),
            ([[x0, y0, z0], [x0, y0, z1], [x1, y0, z1], [x1, y0, z0]], (0, -1, 0)),
            ([[x0, y1, z0], [x1, y1, z0], [x1, y1, z1], [x0, y1, z1]], (0, 1, 0)),
            ([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0]], (0, 0, -1)),
            ([[x0, y0, z1], [x0, y1, z1], [x1, y1, z1], [x1, y0, z1]], (0, 0, 1)),
        )
        factor = scale_factor[label]
        centre = (lower[k] + upper[k]) / 2 if factor != 1 else None
        for face, offset in candidates:
            if aligned:
                neighbour = (label, int(index[k][0]) + offset[0],
                             int(index[k][1]) + offset[1], int(index[k][2]) + offset[2])
                if neighbour in occupied:
                    continue
            if centre is not None:
                face = (centre + (np.asarray(face) - centre) * factor).tolist()
            faces.append(face)
            face_labels.append(label)
    return np.asarray(faces, dtype=float), np.asarray(face_labels, dtype=int), aligned



def _shaded_facecolors(faces, face_labels, cmap, cmap_norm, light_azdeg=300.0,
                       light_altdeg=55.0, strength=0.32, highlight_strength=0.12):
    """Directional shading that keeps each Morse set's colour recognisable.

       Axis-aligned faces otherwise render as flat blocks in which the three
       visible sides of a cube are indistinguishable, so the shape reads as a
       silhouette. Diffuse light separates the sides while the shading stays
       weak enough that the set's colour still identifies it.
    """
    from matplotlib.colors import LightSource, to_rgb

    base = np.asarray([to_rgb(matplotlib.colors.to_hex(cmap(cmap_norm(int(l)))))
                       for l in face_labels], dtype=float)
    if strength == 0:
        return base
    # Face vertices wind inward, so flip the cross-product normals before
    # applying the light.
    normals = -np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = np.divide(normals, lengths, out=np.zeros_like(normals), where=lengths > 0)
    direction = LightSource(azdeg=light_azdeg, altdeg=light_altdeg).direction
    diffuse = np.clip(normals @ direction, 0.0, 1.0)
    shaded = base * (1.0 - strength * (1.0 - diffuse))[:, None]
    highlight = highlight_strength * np.clip((diffuse - 0.4) / 0.6, 0.0, 1.0)
    shaded += (1.0 - shaded) * highlight[:, None]
    return np.clip(shaded, 0.0, 1.0)


def PlotMorseSets3D(morse_sets, morse_nodes=None, cmap=None, clist=None, scale_factor=None,
                    fig_w=8, fig_h=8, xlim=None, ylim=None, zlim=None, axis_labels=True,
                    xlabel='$x$', ylabel='$y$', zlabel='$z$', fontsize=15, elev=22, azim=-55,
                    alpha=1.0, edge_clr=None, linewidth=0.0, lighting=True,
                    light_azdeg=300.0, light_altdeg=55.0, shade_strength=0.32,
                    grid=False, max_ticks=4, zlabel_pos=None, fig_fname=None,
                    dpi=None, rasterize=None, show=None):
    """Plot three-dimensional Morse sets as cubical surfaces.

       Only faces on the boundary of a Morse set are drawn, so the cost follows
       the surface of the sets rather than their volume and a dense block
       renders as a shell. See _exposed_faces for when culling applies.

       lighting shades the faces by their orientation, which is on by default:
       without it the three visible sides of every cube take the same colour
       and the geometry reads as a flat silhouette. Set it False for flat fill.

       grid keeps Matplotlib's grey background panes and their grid lines. Off
       by default: the sets then sit on a plain white ground with only the
       axis lines and ticks, so the surfaces carry the depth cues.

       max_ticks caps the major ticks per axis. A 3-D projection crowds tick
       numbers along three foreshortened edges, so far fewer fit legibly than
       on a plane; None keeps Matplotlib's default density.

       zlabel_pos places the z label in axes fractions. Matplotlib clips a 3-D
       z label when the figure is saved with a tight bounding box, so the label
       is drawn as an unclipped 2-D annotation instead of on the axis.

       rasterize draws the faces as one bitmap inside the vector figure, at
       dpi (600 by default when rasterizing); axes, ticks and labels stay
       vector. None, the default, decides by size. A cubical surface is a
       polygon per exposed cell face, painted in depth order, so unlike the
       2-D plot it cannot be merged into fewer paths: vector output grows
       with the face count while the bitmap grows far more slowly, so sets
       with more than RASTERIZE_FACES faces are rasterized and smaller ones
       stay fully vector. True or False forces either.

       Returns (fig, ax).
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    rows, num_morse_sets = _load_morse_sets(morse_sets)
    dim = _box_dim(rows)
    assert dim == 3, f"PlotMorseSets3D needs 3-D Morse boxes; got dim={dim}"
    num_morse_sets, morse_nodes, scale_factor, cmap, cmap_norm = _resolve_plot_setup(
        rows, num_morse_sets, morse_nodes, cmap, clist, scale_factor)
    morse_nodes = list(morse_nodes)

    faces, face_labels, aligned = _exposed_faces(rows, morse_nodes, scale_factor)
    if len(faces) == 0:
        raise ValueError("no Morse set boxes to plot")
    if rasterize == None:
        rasterize = len(faces) > RASTERIZE_FACES
    if lighting:
        facecolors = _shaded_facecolors(faces, face_labels, cmap, cmap_norm,
                                        light_azdeg=light_azdeg, light_altdeg=light_altdeg,
                                        strength=shade_strength)
    else:
        facecolors = [matplotlib.colors.to_hex(cmap(cmap_norm(int(l))), keep_alpha=True)
                      for l in face_labels]

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    collection = Poly3DCollection(
        faces, facecolors=facecolors,
        edgecolors=('none' if edge_clr == None else edge_clr),
        linewidths=linewidth, alpha=alpha, rasterized=rasterize, shade=False)
    ax.add_collection3d(collection)

    values = np.asarray([[float(v) for v in rect] for rect in rows], dtype=float)
    keep = np.isin(values[:, 6].astype(int), morse_nodes)
    lower, upper = values[keep, :3].min(axis=0), values[keep, 3:6].max(axis=0)
    span = np.maximum(upper - lower, 1e-12)
    margin = 0.04 * span
    ax.set_xlim(list(xlim) if xlim != None else [lower[0] - margin[0], upper[0] + margin[0]])
    ax.set_ylim(list(ylim) if ylim != None else [lower[1] - margin[1], upper[1] + margin[1]])
    ax.set_zlim(list(zlim) if zlim != None else [lower[2] - margin[2], upper[2] + margin[2]])
    ax.view_init(elev=elev, azim=azim)
    if not grid:
        # The grid lines live on the background panes, so hiding one without
        # the other leaves either floating lines or blank grey walls.
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_visible(False)
    if max_ticks != None:
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_major_locator(matplotlib.ticker.MaxNLocator(max_ticks))
    if axis_labels:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        # The projected z label is clipped by a tight bounding box even though
        # it is inside the figure, so keep it off the axis and draw an
        # unclipped copy clear of the tick numbers.
        ax.set_zlabel('')
        x_pos, y_pos = (1.04, 0.55) if zlabel_pos == None else zlabel_pos
        ax.text2D(x_pos, y_pos, zlabel, transform=ax.transAxes, rotation=90,
                  rotation_mode='anchor', ha='center', va='center',
                  fontsize=fontsize, clip_on=False)
    ax.tick_params(labelsize=fontsize)
    return _finish(fig, ax, fig_fname, dpi, show, rasterize)


def PlotMorseSetsScatter(morse_sets, morse_nodes=None, proj_dims=None, cmap=None, clist=None,
                         scale_factor=None, fig_w=8, fig_h=8, xlim=None, ylim=None, margin=0.02,
                         axis_labels=True, xlabel='$x$', ylabel='$y$', fontsize=15,
                         fig_fname=None, dpi=300, rasterize=False, show=None):
    """Plot Morse sets as square scatter markers sized in data units.

       The original Morse set plot, kept under its own name now that
       PlotMorseSets draws rectangles. Markers are square regardless of a box's
       aspect ratio and are sized by the larger side, so they can overlap; that
       is often what you want when boxes are near the resolution limit and
       exact rectangles would vanish.

       Returns (fig, ax).
    """
    rows, num_morse_sets = _load_morse_sets(morse_sets)
    dim = _box_dim(rows)
    if dim == 1:
        # Add extra fake dimension to plot
        rows = [[x0_min, 0, x0_max, x0_max - x0_min, node] for x0_min, x0_max, node in rows]
        dim = 2
    return PlotBoxesScatter(rows, num_morse_sets=num_morse_sets, morse_nodes=morse_nodes,
                            proj_dims=proj_dims, cmap=cmap, clist=clist,
                            scale_factor=scale_factor, fig_w=fig_w, fig_h=fig_h,
                            xlim=xlim, ylim=ylim, axis_labels=axis_labels, xlabel=xlabel,
                            ylabel=ylabel, fontsize=fontsize, fig_fname=fig_fname, dpi=dpi,
                            rasterize=rasterize, show=show)


def PlotBoxesScatter(morse_sets, num_morse_sets=None, morse_nodes=None, proj_dims=None, cmap=None,
                     clist=None, scale_factor=None, fig_w=8, fig_h=8, xlim=None, ylim=None,
                     axis_labels=True, xlabel='$x$', ylabel='$y$', fontsize=15, fig_fname=None,
                     dpi=300, rasterize=False, show=None):
    """Scatter plot of labelled boxes, markers sized in data units."""
    rows = list(morse_sets)
    dim = _box_dim(rows)
    num_morse_sets, morse_nodes, scale_factor, cmap, cmap_norm = _resolve_plot_setup(
        rows, num_morse_sets, morse_nodes, cmap, clist, scale_factor)
    if proj_dims == None:
        d1, d2 = 0, 1
    else:
        d1, d2 = proj_dims[0], proj_dims[1]
    assert max(d1, d2) < dim, "Wrong projection dimensions"

    x_min, x_max, y_min, y_max = _axis_limits(rows, morse_nodes, dim, d1, d2, xlim, ylim,
                                              margin=margin)
    x_axis_width = x_max - x_min
    y_axis_height = y_max - y_min
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    # The scatter plot in Matplotlib uses a marker size in points**2 units.
    # The relationship between points and pixels is 1 point = dpi/72 pixels,
    # hence multiplying by 72 / fig.dpi gives the number of pixels in the plot.
    s0_x = (ax.get_window_extent().width / x_axis_width) * (72.0 / fig.dpi)
    s0_y = (ax.get_window_extent().height / y_axis_height) * (72.0 / fig.dpi)
    for morse_node in morse_nodes:
        morse_set = [rect for rect in rows if int(rect[-1]) == morse_node]
        if not morse_set:
            continue
        clr = matplotlib.colors.to_hex(cmap(cmap_norm(morse_node)), keep_alpha=True)
        X = []; Y = []; S = []
        for rect in morse_set:
            p1 = [rect[d1], rect[d2]]                   # Lower point
            p2 = [rect[dim + d1], rect[dim + d2]]       # Upper point
            p = list((np.array(p1) + np.array(p2)) / 2) # Center point
            s = list((np.array(p2) - np.array(p1)))     # Rect size
            s_x = (scale_factor[morse_node] * s0_x * s[0]) ** 2
            s_y = (scale_factor[morse_node] * s0_y * s[1]) ** 2
            X.append(p[0])
            Y.append(p[1])
            S.append(max(s_x, s_y))
        ax.scatter(X, Y, s=S, marker='s', c=clr, rasterized=rasterize)
    if axis_labels:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    return _finish(fig, ax, fig_fname, dpi, show, rasterize)
