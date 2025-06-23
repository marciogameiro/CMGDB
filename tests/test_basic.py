import CMGDB

def f(x):
    return [x[0] / (2.0 - x[0]), x[1] / (2.0 - x[1])]

def F(rect):
    # Define box map for f
    return CMGDB.BoxMap(f, rect, padding=False)

def test_morse_graph():
    subdiv_min = 6
    subdiv_max = 10
    subdiv_init = 4
    subdiv_limit = 10000
    lower_bounds = [0, 0]
    upper_bounds = [1.2, 1.2]
    model = CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, subdiv_limit, lower_bounds, upper_bounds, F)
    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)
    assert morse_graph.num_vertices() == 4

def test_conley_index():
    subdiv_min = 6
    subdiv_max = 10
    subdiv_init = 4
    subdiv_limit = 10000
    lower_bounds = [0, 0]
    upper_bounds = [1.2, 1.2]
    model = CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, subdiv_limit, lower_bounds, upper_bounds, F)
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(model)
    assert morse_graph.annotations(0) == ['x-1', '0', '0']
    assert morse_graph.annotations(1) == ['0', 'x-1', '0']
    assert morse_graph.annotations(2) == ['0', 'x-1', '0']
    assert morse_graph.annotations(3) == ['0', '0', 'x-1']
