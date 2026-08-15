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

# Vectorized version of BoxMap for use with Model.set_batch_map
def BoxMapBatch(f, rects, mode='corners', padding=False):
    """Evaluate a box map on many rectangles at once.

    f must map an (m, dim) NumPy array of points to an (m, dim) array of
    image points. rects is an (N, 2*dim) array of rectangles, each row
    holding dim lower bounds followed by dim upper bounds. Returns an
    (N, 2*dim) array of image rectangles, equivalent to applying
    BoxMap(f_scalar, rect, mode, padding) to each row.

    Typical use:
        def F_batch(rects):
            return CMGDB.BoxMapBatch(f, rects)
        model = CMGDB.Model(subdiv_min, subdiv_max, lower_bounds, upper_bounds, F)
        model.set_batch_map(F_batch)
    """
    rects = np.asarray(rects, dtype=float)
    N, two_dim = rects.shape
    dim = two_dim // 2
    lower = rects[:, :dim]
    upper = rects[:, dim:]
    if mode == 'corners': # Compute at corner points
        # Evaluate f on all 2^dim corners of every rectangle in one call,
        # then take the componentwise min/max over the corners of each box
        num_corners = 2 ** dim
        corners = np.empty((num_corners, N, dim))
        for k in range(num_corners):
            mask = np.array([(k >> d) & 1 for d in range(dim)], dtype=bool)
            corners[k] = np.where(mask, upper, lower)
        Y = np.asarray(f(corners.reshape(num_corners * N, dim)))
        Y = Y.reshape(num_corners, N, dim)
        Y_lower = Y.min(axis=0)
        Y_upper = Y.max(axis=0)
    elif mode == 'center': # Compute at center point
        padding = True # Must be true for this case
        Y = np.asarray(f((lower + upper) / 2))
        Y_lower = Y
        Y_upper = Y.copy()
    else: # Unknown mode
        raise ValueError("BoxMapBatch supports modes 'corners' and 'center'")
    if padding:
        pad = upper - lower
        Y_lower = Y_lower - pad
        Y_upper = Y_upper + pad
    return np.hstack([Y_lower, Y_upper])

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
        upper_bounds = rect[dim:]
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
