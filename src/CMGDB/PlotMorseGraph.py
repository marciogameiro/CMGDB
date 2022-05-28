# PlotMorseGraph.py
# Marcio Gameiro
# MIT LICENSE
# 2022-05-28

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import graphviz

def PlotMorseGraph(morse_graph, cmap=None, clist=None, shape=None, margin=None):
    """Plot Morse graph using cmap as the colormap or clist
       as a list of colors to make a colormap."""

    # Default colormap
    default_cmap = matplotlib.cm.brg
    # Number of vertices
    num_verts = len(morse_graph.vertices())
    # Set defaults for unset values
    if clist and cmap == None:
        cmap = matplotlib.colors.ListedColormap(clist[:num_verts], name='clist')
    if cmap == None:
        cmap = default_cmap
    if shape == None:
        shape = 'ellipse'
    if margin == None:
        margin = '0.11, 0.055'
    margin = str(margin)

    # Normalization for color map
    cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=num_verts-1)

    # Get list of attractors (nodes without children)
    attractors = [v for v in morse_graph.vertices() if len(morse_graph.adjacencies(v)) == 0]

    def graphviz_string(graph):
        # Return graphviz string describing the graph

        def vertex_color(v):
            # Return vertex color
            clr = matplotlib.colors.to_hex(cmap(cmap_norm(v)), keep_alpha=True)
            return str(clr)

        def vertex_label(v):
            """Return vertex label"""
            # Handle the case of annotations not present
            annotations = getattr(graph, 'annotations', [])
            if annotations and len(annotations(v)):
                return str(v) + ' : (' + ', '.join(annotations(v)) + ')'
            return str(v) # Just use index if no annotations

        # Make graphviz string
        gv = 'digraph {\n'
        for v in graph.vertices():
            gv += str(v) + ' [label="' + vertex_label(v) + '"' + (
                ', shape=' + shape + ', style=filled, fillcolor="' + vertex_color(v) +
                '", margin="' + margin + '"];\n')

        # Set rank for attractors
        gv += '{rank=same; '
        for v in attractors:
            gv += str(v) + ' '
        gv += '}; \n'

        # Set the graph edges
        for u in graph.vertices():
            for v in graph.adjacencies(u):
                gv += str(u) + ' -> ' + str(v) + ';\n'

        gv += '}\n' # Close bracket
        return gv

    gv = graphviz_string(morse_graph) # Get graphviz string
    return graphviz.Source(gv) # Return graphviz render
