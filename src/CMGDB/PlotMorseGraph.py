import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import graphviz

def PlotMorseGraph(morse_graph, cmap=matplotlib.cm.brg):
    # Plot Morse graph using cmap as the colormap

    num_verts = len(morse_graph.vertices()) # Number of vertices
    # Normalization for color map
    cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=num_verts-1)

    # Get list of attractors (nodes without children)
    attractors = [v for v in morse_graph.vertices() if len(morse_graph.adjacencies(v)) == 0]

    def graphviz_string(graph):
        # Return graphviz string describing the graph

        def vertex_color(v):
            # Return vertex color
            clr = matplotlib.colors.to_hex(cmap(cmap_norm(v)))
            return str(clr)

        def vertex_label(v):
            # Return vertex label
            if len(graph.annotations(v)):
                return str(v) + ' : (' + ', '.join(graph.annotations(v)) + ')'
            return str(v) # Just use index if no annotations

        # Make graphviz string
        gv = 'digraph {\n'
        for v in graph.vertices():
            gv += str(v) + ' [label="' + vertex_label(v) + '"' + (
                ', style=filled, fillcolor="' + vertex_color(v) + '"];\n')

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
