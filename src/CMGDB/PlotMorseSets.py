import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import CMGDB

def PlotMorseSets(morse_sets, morse_nodes=None, proj_dims=None,
                  cmap=matplotlib.cm.brg, fig_w=8, fig_h=8):
    # Check if morse_sets is a Morse graph, file name, or list
    if type(morse_sets) == CMGDB._cmgdb.MorseGraph: # Morse graph
        morse_graph = morse_sets
        num_morse_sets = morse_graph.num_vertices()
        if morse_nodes == None:
            morse_nodes = range(num_morse_sets)
        morse_sets = [box + [node] for node in morse_nodes for box in morse_graph.morse_set_boxes(node)]
    elif type(morse_sets) == str: # String representing file name
        morse_fname = morse_sets
        morse_sets = CMGDB.LoadMorseSetFile(morse_fname)
        num_morse_sets = None
    else: # List of Morse sets
        num_morse_sets = None
    # Plot Morse set boxes as a scatter plot
    PlotBoxesScatter(morse_sets, num_morse_sets=num_morse_sets, morse_nodes=morse_nodes,
                     proj_dims=proj_dims, cmap=cmap, fig_w=fig_w, fig_h=fig_h)

def PlotBoxesScatter(morse_sets, num_morse_sets=None, morse_nodes=None,
                     proj_dims=None, cmap=matplotlib.cm.brg, fig_w=8, fig_h=8):
    rect = morse_sets[0]
    assert len(rect) % 2 == 1, "Wrong dimension in Morse sets data"
    dim = int((len(rect) - 1) / 2)
    if num_morse_sets == None:
        num_morse_sets = max([int(rect[-1]) for rect in morse_sets]) + 1
    if morse_nodes == None:
        morse_nodes = range(num_morse_sets)
    if proj_dims == None:
        d1 = 0
        d2 = 1
    else:
        d1 = proj_dims[0]
        d2 = proj_dims[1]
    assert max(d1, d2) < dim, "Wrong projection dimensions"
    # Get min and max x values
    x_lower = [rect[d1] for rect in morse_sets if int(rect[-1]) in morse_nodes]
    x_upper = [rect[dim + d1] for rect in morse_sets if int(rect[-1]) in morse_nodes]
    x_axis_width = max(x_upper) - min(x_lower)
    # Next plot Morse sets
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
    # Normalization for color map
    cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=num_morse_sets-1)
    # Get unit size in data units for marker size
    s0 = ((ax.get_window_extent().width / x_axis_width) * (72.0 / fig.dpi)) ** 2
    for morse_node in morse_nodes:
        morse_set = [rect for rect in morse_sets if rect[-1] == morse_node]
        # Use morse_node as color index for consistency if not plotting all
        clr = matplotlib.colors.to_hex(cmap(cmap_norm(morse_node)))
        X = []; Y = []; S = []
        for rect in morse_set:
            p1 = [rect[d1], rect[d2]]                   # Lower point
            p2 = [rect[dim + d1], rect[dim + d2]]       # Upper point
            p = list((np.array(p1) + np.array(p2)) / 2) # Center point
            s = list((np.array(p2) - np.array(p1)))     # Rect size
            X.append(p[0])
            Y.append(p[1])
            S.append(s0 * max(s) ** 2)
        plt.scatter(X, Y, s=S, marker='s', c=clr)
        # plt.scatter(X, Y, s=S, marker='s', c=clr, alpha=0.5)
    plt.show()
