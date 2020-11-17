import numpy as np
import csv
import CMGDB

def SaveMorseSets(morse_graph, morse_graph_fname):
    num_morse_sets = morse_graph.num_vertices()
    with open(morse_graph_fname, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        for node in range(num_morse_sets):
            for morse_box in morse_graph.morse_set_boxes(node):
                morse_rect = morse_box + [node]
                csv_writer.writerow(morse_rect)
