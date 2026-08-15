### BoxMapData.py
### MIT LICENSE 2024 Marcio Gameiro

import numpy as np

class BoxMapData:
    """Define a box map from datasets X and Y, where the points in Y are the images
       of the points in X. Given an input rectangle rect, compute a rectangle f_rect
       containing the points in Y which are image of the points in X inside of rect.
       If there are no X points in rect the result depends on the flag map_empty: If
       map_empty is 'outside' return a rectangle outiside the domain, if map_empty is
       'terminate' raise an exception, and if map_empty is 'interp' use a form of
       interpolation to compute the image.

       The points inside a rectangle are found with a uniform-grid spatial index
       built once over X: each query gathers candidate points from only the grid
       bins the rectangle overlaps, then applies the exact containment test to the
       candidates. The selected points -- and therefore all results -- are identical
       to the linear-scan reference implementation BoxMapDataLinear, but each query
       costs O(candidates) instead of a scan of the full dataset.

       With use_index='auto' (the default) small datasets skip the index and scan
       directly, since below about a thousand points a single vectorized scan is
       cheaper than the index bookkeeping. Pass use_index=True/False to force
       either behavior; results are identical in all cases."""

    _INDEX_THRESHOLD = 1000

    def __init__(self, X, Y, map_empty='interp', lower_bounds=None, upper_bounds=None,
                 domain_padding=True, padding=False, points_per_bin=8, use_index='auto'):
        if map_empty not in ['interp', 'outside', 'terminate']:
            raise ValueError("Invalid value for map_empty. Allowed values are: 'interp', 'outside', or 'terminate'")
        if map_empty == 'outside' and (lower_bounds is None or upper_bounds is None):
            raise ValueError("The bounds lower_bounds and upper_bounds must be provided if map_empty is 'outside'")
        self.X = np.array(X)
        self.Y = np.array(Y)
        self.map_empty = map_empty
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.domain_padding = domain_padding
        self.padding = padding
        # num_pts, dim = self.X.shape
        self.dim = self.X.shape[1]
        if use_index == 'auto':
            use_index = self.X.shape[0] > self._INDEX_THRESHOLD
        self._use_index = bool(use_index)
        if self._use_index:
            self._build_index(points_per_bin)

    def _build_index(self, points_per_bin):
        """Build the uniform-grid index over X: assign each point an integer bin,
           sort point indices by bin, and store CSR-style bin offsets so the points
           of bin b are self._sorted_idx[self._bin_offsets[b]:self._bin_offsets[b+1]].
           The number of bins targets points_per_bin points per bin on average."""
        num_pts = self.X.shape[0]
        self._bins_per_dim = max(1, int(round((num_pts / points_per_bin) ** (1.0 / self.dim))))
        self._x_min = self.X.min(axis=0).astype(float)
        x_max = self.X.max(axis=0).astype(float)
        # Guard against zero-width dimensions (all points equal in a coordinate)
        widths = x_max - self._x_min
        widths[widths == 0] = 1.0
        self._bin_widths = widths / self._bins_per_dim
        # Integer bin coordinates per point, clamped so boundary points fall in the last bin
        coords = ((self.X - self._x_min) / self._bin_widths).astype(np.int64)
        np.clip(coords, 0, self._bins_per_dim - 1, out=coords)
        # Flatten multi-dimensional bin coordinates to a single bin id
        self._bin_multipliers = self._bins_per_dim ** np.arange(self.dim, dtype=np.int64)
        bin_ids = coords @ self._bin_multipliers
        self._num_bins = self._bins_per_dim ** self.dim
        self._sorted_idx = np.argsort(bin_ids, kind='stable')
        sorted_bins = bin_ids[self._sorted_idx]
        self._bin_offsets = np.searchsorted(sorted_bins, np.arange(self._num_bins + 1))

    def _candidate_indices(self, l_bounds, u_bounds):
        """Return the (sorted) indices of all points in bins overlapping the given
           rectangle. This is a superset of the points inside the rectangle; the
           exact containment test is applied by the caller. Rectangles reaching
           outside the data range clamp to the edge bins."""
        lo = np.floor((np.asarray(l_bounds, dtype=float) - self._x_min) / self._bin_widths).astype(np.int64)
        hi = np.floor((np.asarray(u_bounds, dtype=float) - self._x_min) / self._bin_widths).astype(np.int64)
        np.clip(lo, 0, self._bins_per_dim - 1, out=lo)
        np.clip(hi, 0, self._bins_per_dim - 1, out=hi)
        # Rectangles covering every bin (e.g. from interpolate's doubling loop)
        # degrade to the full dataset
        if np.all(lo == 0) and np.all(hi == self._bins_per_dim - 1):
            return np.arange(self.X.shape[0])
        # All bin ids in the hyper-rectangle of bins [lo, hi]
        axes = [np.arange(lo[d], hi[d] + 1) * self._bin_multipliers[d] for d in range(self.dim)]
        bin_ids = axes[0]
        for axis in axes[1:]:
            bin_ids = (bin_ids[:, None] + axis[None, :]).ravel()
        parts = [self._sorted_idx[self._bin_offsets[b]:self._bin_offsets[b + 1]] for b in bin_ids]
        candidates = np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)
        # Sort so downstream results keep the dataset's row order, making the
        # output identical to the linear-scan boolean mask
        return np.sort(candidates)

    def __call__(self, rect):
        return self.compute(rect)

    def map_points(self, rect):
        """Return the points in Y which are image of the points in X inside of rect."""
        l_bounds = rect[:self.dim]
        u_bounds = rect[self.dim:]
        if not self._use_index:
            # Small dataset: a single vectorized scan is cheaper than the index
            index_mask = np.all((self.X >= l_bounds) & (self.X <= u_bounds), axis=1)
            return self.Y[index_mask]
        # Candidate points from the bins overlapping rect
        candidates = self._candidate_indices(l_bounds, u_bounds)
        X_cand = self.X[candidates]
        # Exact containment test on the candidates (same test as the linear scan)
        index_mask = np.all((X_cand >= l_bounds) & (X_cand <= u_bounds), axis=1)
        # Get the corresponding points in Y
        Y_rect = self.Y[candidates[index_mask]]
        return Y_rect

    def interpolate(self, rect):
        """Compute the image of the empty rectangle rect using interpolation.
           Double the size of the rectangle until there are X points inside
           and return a rectangle with the corresponding points in Y."""
        Y_rect = self.map_points(rect)
        if Y_rect.size > 0:
            return Y_rect
        l_bounds = rect[:self.dim]
        u_bounds = rect[self.dim:]
        # Double rectangle size until nonempty
        while Y_rect.size == 0:
            l_bounds_new = [l_b - (u_b - l_b) / 2 for l_b, u_b in zip(l_bounds, u_bounds)]
            u_bounds_new = [u_b + (u_b - l_b) / 2 for l_b, u_b in zip(l_bounds, u_bounds)]
            l_bounds, u_bounds = l_bounds_new, u_bounds_new
            rect_new = l_bounds_new + u_bounds_new
            Y_rect = self.map_points(rect_new)
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
        # Pad domain rectangle if flag is set
        if self.domain_padding:
            l_bounds = [l_b - (u_b - l_b) for l_b, u_b in zip(rect[:self.dim], rect[self.dim:])]
            u_bounds = [u_b + (u_b - l_b) for l_b, u_b in zip(rect[:self.dim], rect[self.dim:])]
            rect_new = l_bounds + u_bounds
            Y_rect = self.map_points(rect_new)
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

    def batch(self, rects):
        """Evaluate the box map on many rectangles at once. Suitable for use with
           Model.set_batch_map: model.set_batch_map(F.batch)."""
        return [self.compute(list(rect)) for rect in np.asarray(rects, dtype=float)]
