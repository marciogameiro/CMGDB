### BoxMapData.py
### MIT LICENSE 2024 Marcio Gameiro

# TODO: Use a KD-Tree or a grid to make it more efficient to find the points inside a rectangle

import numpy as np

class BoxMapData:
    """Define a box map from datasets X and Y, where the points in Y are the images
       of the points in X. Given an input rectangle rect, compute a rectangle f_rect
       containing the points in Y which are image of the points in X inside of rect.
       If there are no X points in rect the result depends on the flag map_empty: If
       map_empty is 'outside' return a rectangle outiside the domain, if map_empty is
       'terminate' raise an exception, and if map_empty is 'interp' use a form of
       interpolation to compute the image."""

    def __init__(self, X, Y, map_empty='interp', lower_bounds=None, upper_bounds=None, padding=False):
        if map_empty not in ['interp', 'outside', 'terminate']:
            raise ValueError("Invalid value for map_empty. Allowed values are: 'interp', 'outside', or 'terminate'")
        if map_empty == 'outside' and (lower_bounds is None or upper_bounds is None):
            raise ValueError("The bounds lower_bounds and upper_bounds must be provided if map_empty is 'outside'")
        self.X = np.array(X)
        self.Y = np.array(Y)
        self.map_empty = map_empty
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.padding = padding
        # num_pts, dim = self.X.shape
        self.dim = self.X.shape[1]

    def __call__(self, rect):
        return self.compute(rect)

    def map_points(self, rect):
        """Return the points in Y which are image of the points in X inside of rect."""
        l_bounds = rect[:self.dim]
        u_bounds = rect[self.dim:]
        # Get index mask for the points in X inside rect
        index_mask = np.all((self.X >= l_bounds) & (self.X <= u_bounds), axis=1)
        # Get the corresponding points in Y
        Y_rect = self.Y[index_mask]
        return Y_rect

    def interpolate(self, rect):
        """Compute the image of the empty rectangle rect using interpolation. Double
           the size of the rectangle until there are X points inside, then return a
           rectangle with the corresponding points in Y (scale down the size of the
           rectangle proportionally to the increase in size of the input rectangle)."""
        Y_rect = self.map_points(rect)
        if Y_rect.size > 0:
            return Y_rect
        l_bounds = rect[:self.dim]
        u_bounds = rect[self.dim:]
        n_steps = 0
        # Double rectangle size until nonempty
        while Y_rect.size == 0:
            l_bounds_new = [l_b - (u_b - l_b) / 2 for l_b, u_b in zip(l_bounds, u_bounds)]
            u_bounds_new = [u_b + (u_b - l_b) / 2 for l_b, u_b in zip(l_bounds, u_bounds)]
            l_bounds, u_bounds = l_bounds_new, u_bounds_new
            rect_new = l_bounds_new + u_bounds_new
            Y_rect = self.map_points(rect_new)
            n_steps += 1
        # Factor by which rect was increased
        factor = 2**n_steps
        # Get min and max vertices of Y_rect
        Y_l_bounds = list(np.min(Y_rect, axis=0))
        Y_u_bounds = list(np.max(Y_rect, axis=0))
        # Get center and radius (half width) of Y_rect
        Y_center_radius = [((l_b + u_b) / 2, (u_b - l_b) / 2) for l_b, u_b in zip(Y_l_bounds, Y_u_bounds)]
        # Return min and max vertices of scaled-down rectangle
        y_l_bounds = [c - r / factor for c, r in Y_center_radius]
        y_u_bounds = [c + r / factor for c, r in Y_center_radius]
        Y_rect = np.array([y_l_bounds, y_u_bounds])
        return Y_rect

    def compute(self, rect):
        """Compute a rectangle containing the images of the points inside rect."""
        # Compute points in the image
        Y_rect = self.map_points(rect)
        # Raise an exception if empty image and map_empty is terminate
        if Y_rect.size == 0 and self.map_empty == 'terminate':
            raise ValueError(f'Rectangle {rect} has empty image')
        # Map to a box outside if empty image and map_empty is outside
        if Y_rect.size == 0 and self.map_empty == 'outside':
            # Make a box outside the domain
            f_l_bounds = [b + 1 for b in self.upper_bounds]
            f_u_bounds = [b + 2 for b in self.upper_bounds]
            f_rect = f_l_bounds + f_u_bounds
            return f_rect
        # Interpolate if empty image and map_empty is interp
        if Y_rect.size == 0 and self.map_empty == 'interp':
            Y_rect = self.interpolate(rect)
        # Get lower and upper bounds of Y_rect
        Y_l_bounds = list(np.min(Y_rect, axis=0))
        Y_u_bounds = list(np.max(Y_rect, axis=0))
        # Add padding if padding is True
        f_l_bounds = [Y_l_bounds[k] - (rect[k + self.dim] - rect[k] if self.padding else 0) for k in range(self.dim)]
        f_u_bounds = [Y_u_bounds[k] + (rect[k + self.dim] - rect[k] if self.padding else 0) for k in range(self.dim)]
        f_rect = f_l_bounds + f_u_bounds
        return f_rect
