import numpy as np
import itertools
import CMGDB

# Returns corner points of a rectangle
def CornerPoints(rect):
    dim = int(len(rect) / 2)
    # Get list of intervals
    list_intvals = [[rect[d], rect[d + dim]] for d in range(dim)]
    # Get points in the cartesian product of intervals
    X = [list(u) for u in itertools.product(*list_intvals)]
    return X

# Returns center point of a rectangle
def CenterPoint(rect):
    dim = int(len(rect) / 2)
    x_center = [(rect[d] + rect[dim + d]) / 2 for d in range(dim)]
    return [x_center]

# Return sample points in rectangle
def SamplePoints(lower_bounds, upper_bounds, num_pts):
    # Sample num_pts in dimension dim, where each
    # component of the sampled points are in the
    # ranges given by lower_bounds and upper_bounds
    dim = len(lower_bounds)
    X = np.random.uniform(lower_bounds, upper_bounds, size=(num_pts,dim))
    return list(X)

# Map that takes a rectangle and returns a rectangle
def BoxMap(f, rect, mode='corners', padding=False, num_pts=10):
    dim = int(len(rect) / 2)
    if mode == 'corners': # Compute at corner points
        X = CornerPoints(rect)
    elif mode == 'center': # Compute at center point
        padding = True # Must be true for this case
        X = CenterPoint(rect)
    elif mode == 'random': # Compute at random point
        # Get lower and upper bounds
        lower_bounds = rect[:dim]
        upper_nounds = rect[dim:]
        X = SamplePoints(lower_bounds, upper_bounds, num_pts)
    else: # Unknown mode
        return []
    # Evaluate f at point in X
    Y = [f(x) for x in X]
    # Get lower and upper bounds of Y
    Y_l_bounds = [min([y[d] for y in Y]) - ((rect[d + dim] - rect[d]) if padding else 0) for d in range(dim)]
    Y_u_bounds = [max([y[d] for y in Y]) + ((rect[d + dim] - rect[d]) if padding else 0) for d in range(dim)]
    f_rect = Y_l_bounds + Y_u_bounds
    return f_rect
