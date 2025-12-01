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
    # Default color list
    default_clist = ['#1f77b4', '#e6550d', '#31a354', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                     '#bcbd22', '#80b1d3', '#ffffb3', '#fccde5', '#b3de69', '#fdae6b', '#6a3d9a', '#c49c94',
                     '#fb8072', '#dbdb8d', '#bc80bd', '#ffed6f', '#637939', '#c5b0d5', '#636363', '#c7c7c7',
                     '#8dd3c7', '#b15928', '#e8cb32', '#9e9ac8', '#74c476', '#ff7f0e', '#9edae5', '#90d743',
                     '#e7969c', '#17becf', '#7b4173', '#8ca252', '#ad494a', '#8c6d31', '#a55194', '#00cc49']
    # # Default colormap
    # default_cmap = matplotlib.cm.tab20
    # Number of vertices
    num_verts = len(morse_graph.vertices())
    # Set defaults for unset values
    if shape == None:
        shape = 'ellipse'
    if margin == None:
        margin = '0.11, 0.055'
    margin = str(margin)
    # Set colormap for Morse graph
    if cmap == None and clist == None:
        clist = default_clist
    if cmap == None:
        cmap = matplotlib.colors.ListedColormap(clist[:num_verts])
    # Get number of colors in the colormap
    try:
        # Colormap is listed colors colormap
        num_colors = len(cmap.colors)
    except:
        num_colors = 0 # Colormap is "continuous"
    if (num_colors > 0) and (num_colors < num_verts):
        # Make colormap cyclic
        cmap_norm = lambda k: k % num_colors
    else:
        # Normalization for color map
        cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=num_verts-1)
    # Get list of attractors (nodes without children)
    attractors = [v for v in morse_graph.vertices() if len(morse_graph.adjacencies(v)) == 0]
    # Get set of children nodes (nodes with a parent)
    children_nodes = set()
    for v in morse_graph.vertices():
        children_nodes.update(morse_graph.adjacencies(v))
    # Get list of repellers (nodes without parents)
    repellers = list(set(morse_graph.vertices()).difference(children_nodes))

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

        # Set rank for repellers
        gv += '{rank=same; '
        for v in repellers:
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
