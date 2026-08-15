### PrecomputedBoxMap.py
### MIT LICENSE 2026 Marcio Gameiro

import sys
import numpy as np

class PrecomputedBoxMap:
    """Box map backed by a table of map evaluations precomputed on the corner
       lattice of the finest subdivision grid.

       The constructor evaluates the map f once at every corner node of the
       dyadic grid at depth subdiv_max, in memory-bounded chunks. After that,
       calling the object with a rectangle performs no map evaluation at all:
       the corners of any box of the subdivision tree (at any depth up to
       subdiv_max) lie exactly on the precomputed lattice, so the image
       rectangle is a table lookup followed by a componentwise min/max --
       the same combinatorial box image BoxMap(f, rect) would produce.

       This pays off when f is expensive (neural network surrogates, Gaussian
       processes, ODE integration): each lattice point is evaluated exactly
       once even though adjacent boxes share it and it recurs across every
       subdivision level, and the evaluation happens in large uniform batches.
       For cheap analytic maps the live BoxMap / BoxMapBatch path is simpler
       and just as fast. Memory scales with the lattice (about 2**subdiv_max
       nodes), so this suits moderate depths.

       f must map an (m, dim) NumPy array of points to an (m, dim) array of
       image points (the BoxMapBatch convention). If PyTorch is in use and f
       is a torch.nn.Module, it is evaluated in float32 on device ('auto'
       selects mps, then cuda, then cpu); torch is never imported otherwise.

       Typical use:
           F = CMGDB.PrecomputedBoxMap(f, lower_bounds, upper_bounds, subdiv_max)
           model = CMGDB.Model(subdiv_min, subdiv_max, lower_bounds, upper_bounds, F)
           model.set_batch_map(F.batch)
    """

    def __init__(self, f, lower_bounds, upper_bounds, subdiv_max,
                 mode='corners', padding=False, batch_points='auto', device='auto'):
        self.lower_bounds = np.asarray(lower_bounds, dtype=float)
        self.upper_bounds = np.asarray(upper_bounds, dtype=float)
        self.dim = self.lower_bounds.shape[0]
        if self.upper_bounds.shape != (self.dim,):
            raise ValueError("lower_bounds and upper_bounds must have the same length")
        if np.any(self.upper_bounds <= self.lower_bounds):
            raise ValueError("upper_bounds must be strictly greater than lower_bounds")
        subdiv_max = int(subdiv_max)
        if subdiv_max < 1:
            raise ValueError(f"subdiv_max must be positive; got {subdiv_max}")
        self.subdiv_max = subdiv_max

        if mode == 'corners':
            # Sample each box at its 2^dim corners, which lie on the plain
            # corner lattice (scale 1)
            scale = 1
            numerators = np.array([[(k >> d) & 1 for d in range(self.dim)]
                                   for k in range(2 ** self.dim)], dtype=np.int64)
        elif mode == 'center':
            # Sample each box at its center, which lies on the lattice refined
            # by one extra level (scale 2). Center mode must pad (as in BoxMap)
            padding = True
            scale = 2
            numerators = np.array([[1] * self.dim], dtype=np.int64)
        else:
            raise ValueError("PrecomputedBoxMap supports modes 'corners' and 'center'")
        self._scale = scale
        self._numerators = numerators
        self.padding = padding

        # CMGDB bisects coordinate (depth % dim) at each depth, so after
        # subdiv_max subdivisions axis j has been split ceil((subdiv_max-j)/dim)
        # times; using the per-axis counts (rather than the max) keeps the
        # table smaller whenever subdiv_max % dim != 0
        axis_depths = [(subdiv_max - j + self.dim - 1) // self.dim for j in range(self.dim)]
        self._cells_per_axis = np.array([2 ** t for t in axis_depths], dtype=np.int64)
        self._finest_box_side = (self.upper_bounds - self.lower_bounds) / self._cells_per_axis
        nodes_per_axis = self._cells_per_axis * scale + 1
        self._nodes_per_axis = nodes_per_axis

        n_total = int(np.prod(nodes_per_axis))
        if n_total > 2 ** 31:
            raise ValueError(
                f"The corner lattice at subdiv_max={subdiv_max} has {n_total} nodes, "
                "which is too large to precompute; reduce subdiv_max or use the "
                "live BoxMap / BoxMapBatch evaluation instead")
        if batch_points == 'auto':
            chunk = min(n_total, 1 << 20)
        else:
            chunk = max(1, int(batch_points))

        evaluator = self._make_evaluator(f, device)
        step = self._finest_box_side / scale
        table = np.empty((n_total, self.dim), dtype=float)
        for start in range(0, n_total, chunk):
            stop = min(start + chunk, n_total)
            flat_idx = np.arange(start, stop, dtype=np.int64)
            multi_idx = np.stack(np.unravel_index(flat_idx, tuple(nodes_per_axis)), axis=-1)
            points = self.lower_bounds + multi_idx * step
            values = np.asarray(evaluator(points), dtype=float)
            if values.shape != (stop - start, self.dim):
                raise ValueError(
                    f"The map must return an array of shape (m, {self.dim}); "
                    f"got {values.shape} for m={stop - start}")
            table[start:stop] = values
        self._table = table.reshape(tuple(nodes_per_axis) + (self.dim,))

    def _make_evaluator(self, f, device):
        # Torch is used only if the caller already imported it and f is a Module
        torch = sys.modules.get('torch')
        if torch is not None and isinstance(f, torch.nn.Module):
            if device == 'auto':
                if torch.backends.mps.is_available():
                    device = 'mps'
                elif torch.cuda.is_available():
                    device = 'cuda'
                else:
                    device = 'cpu'
            module = f.to(device)
            module.eval()

            def evaluator(points):
                with torch.no_grad():
                    x = torch.as_tensor(points, dtype=torch.float32, device=device)
                    return module(x).detach().cpu().numpy().astype(float)
            return evaluator
        return lambda points: f(points)

    def _box_indices(self, lower, upper):
        """Lattice cell indices of one or many boxes. Box bounds produced by the
           subdivision tree are dyadic up to floating point rounding, so rounding
           to the nearest cell index is exact; a large deviation or a zero span
           means the box does not live on the subdiv_max lattice (typically the
           model subdivides deeper than subdiv_max) and is an error."""
        raw_lower = (lower - self.lower_bounds) / self._finest_box_side
        raw_upper = (upper - self.lower_bounds) / self._finest_box_side
        i_lower = np.round(raw_lower).astype(np.int64)
        i_upper = np.round(raw_upper).astype(np.int64)
        deviation = max(np.abs(raw_lower - i_lower).max(), np.abs(raw_upper - i_upper).max())
        if deviation > 0.01:
            raise ValueError(
                "Box bounds do not lie on the subdiv_max lattice "
                f"(deviation {deviation:.3g} cells). Is subdiv_max={self.subdiv_max} "
                "at least the model's maximum subdivision depth?")
        np.clip(i_lower, 0, self._cells_per_axis, out=i_lower)
        np.clip(i_upper, 0, self._cells_per_axis, out=i_upper)
        if np.any(i_upper <= i_lower):
            raise ValueError(
                "Encountered a box finer than the subdiv_max lattice. "
                f"Is subdiv_max={self.subdiv_max} at least the model's "
                "maximum subdivision depth?")
        return i_lower, i_upper

    def __call__(self, rect):
        rect_arr = np.asarray(rect, dtype=float)
        i_lower, i_upper = self._box_indices(rect_arr[:self.dim], rect_arr[self.dim:])
        span = i_upper - i_lower
        # Table node indices of the box's sample points: offset k/scale of the
        # box sits k*span nodes above the box's lower corner (integer at every
        # box depth because the offsets are dyadic)
        nodes = i_lower[None, :] * self._scale + self._numerators * span[None, :]
        samples = self._table[tuple(nodes[:, d] for d in range(self.dim))]
        image_lower = samples.min(axis=0)
        image_upper = samples.max(axis=0)
        if self.padding:
            pad = rect_arr[self.dim:] - rect_arr[:self.dim]
            image_lower = image_lower - pad
            image_upper = image_upper + pad
        return list(image_lower) + list(image_upper)

    def batch(self, rects):
        """Evaluate the box map on many rectangles at once. Suitable for use
           with Model.set_batch_map: model.set_batch_map(F.batch)."""
        rects_arr = np.asarray(rects, dtype=float)
        i_lower, i_upper = self._box_indices(rects_arr[:, :self.dim], rects_arr[:, self.dim:])
        span = i_upper - i_lower
        nodes = (i_lower[:, None, :] * self._scale +
                 self._numerators[None, :, :] * span[:, None, :])
        samples = self._table[tuple(nodes[:, :, d] for d in range(self.dim))]
        image_lower = samples.min(axis=1)
        image_upper = samples.max(axis=1)
        if self.padding:
            pad = rects_arr[:, self.dim:] - rects_arr[:, :self.dim]
            image_lower = image_lower - pad
            image_upper = image_upper + pad
        return np.hstack([image_lower, image_upper])
