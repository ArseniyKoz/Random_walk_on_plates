from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


@dataclass(frozen=True)
class Polyhedron:
    """Convex polyhedron represented as intersection of half-spaces.

    The boundary planes are represented by outward unit normals ``nu`` and
    offsets ``b`` such that each plane is ``dot(nu_i, x) = b_i`` and the
    interior satisfies ``dot(nu_i, x) - b_i <= 0`` for every face.
    """

    nu: FloatArray  # shape (M, 3)
    b: FloatArray  # shape (M,)

    def __post_init__(self) -> None:
        nu = np.asarray(self.nu, dtype=float)
        b = np.asarray(self.b, dtype=float)

        if nu.ndim != 2 or nu.shape[1] != 3:
            raise ValueError("nu must have shape (M, 3).")
        if b.ndim != 1 or b.shape[0] != nu.shape[0]:
            raise ValueError("b must have shape (M,) and match nu.")
        if nu.shape[0] == 0:
            raise ValueError("Polyhedron must contain at least one face.")

        norms = np.linalg.norm(nu, axis=1)
        if np.any(norms <= 0.0):
            raise ValueError("Each face normal must have non-zero norm.")
        if not np.allclose(norms, 1.0, rtol=1e-10, atol=1e-12):
            raise ValueError("All face normals must be unit vectors.")

        object.__setattr__(self, "nu", nu)
        object.__setattr__(self, "b", b)

    @property
    def num_faces(self) -> int:
        return int(self.b.shape[0])

    def signed_distances(self, x: ArrayLike) -> FloatArray:
        """Return signed distances d_i(x)=dot(nu_i,x)-b_i for all faces."""
        x_arr = _as_vec3("x", x)
        return self.nu @ x_arr - self.b

    def is_inside_or_on(self, x: ArrayLike, eps_in: float = 0.0) -> bool:
        """Check whether x satisfies d_i(x) <= eps_in for every face."""
        d = self.signed_distances(x)
        return bool(np.all(d <= eps_in))

    def closest_outside_face_index(self, x: ArrayLike, eps_in: float = 0.0) -> int:
        """Return argmin_{i: d_i(x)>eps_in} d_i(x) for a point outside."""
        d = self.signed_distances(x)
        mask = d > eps_in
        if not np.any(mask):
            raise ValueError("Point is not outside with respect to eps_in.")
        d_masked = np.where(mask, d, np.inf)
        return int(np.argmin(d_masked))

    def characteristic_length(self) -> float:
        """Heuristic geometric scale used for numeric tolerances."""
        scale = float(np.max(np.abs(self.b))) if self.b.size else 1.0
        return max(1.0, scale)


def build_polyhedron_from_planes(
    planes: Sequence[tuple[ArrayLike, ArrayLike]],
) -> Polyhedron:
    """Build ``Polyhedron`` from face list (p_i, nu_i).

    Parameters
    ----------
    planes:
        Sequence of pairs ``(p_i, nu_i)`` where ``p_i`` is any point on the
        i-th face plane and ``nu_i`` is the outward normal.
    """
    if len(planes) == 0:
        raise ValueError("At least one plane is required.")

    nu_list: list[FloatArray] = []
    b_list: list[float] = []
    for idx, (p_i, nu_i) in enumerate(planes):
        p = _as_vec3(f"p[{idx}]", p_i)
        nu = _as_vec3(f"nu[{idx}]", nu_i)
        n_norm = float(np.linalg.norm(nu))
        if n_norm == 0.0:
            raise ValueError(f"nu[{idx}] has zero norm.")
        nu = nu / n_norm
        b = float(np.dot(nu, p))
        nu_list.append(nu)
        b_list.append(b)

    return Polyhedron(nu=np.vstack(nu_list), b=np.asarray(b_list, dtype=float))


def orient_normals(poly: Polyhedron, interior_point: ArrayLike) -> Polyhedron:
    """Orient normals so that interior_point satisfies d_i(c)<=0 for all faces."""
    c = _as_vec3("interior_point", interior_point)
    nu = poly.nu.copy()
    b = poly.b.copy()

    d = nu @ c - b
    flip_mask = d > 0.0
    nu[flip_mask] *= -1.0
    b[flip_mask] *= -1.0
    return Polyhedron(nu=nu, b=b)


def make_axis_aligned_box(min_corner: ArrayLike, max_corner: ArrayLike) -> Polyhedron:
    """Create a box polyhedron represented by six planes only."""
    mn = _as_vec3("min_corner", min_corner)
    mx = _as_vec3("max_corner", max_corner)
    if not np.all(mx > mn):
        raise ValueError("Each component of max_corner must exceed min_corner.")

    planes = [
        (np.array([mn[0], 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])),
        (np.array([mx[0], 0.0, 0.0]), np.array([1.0, 0.0, 0.0])),
        (np.array([0.0, mn[1], 0.0]), np.array([0.0, -1.0, 0.0])),
        (np.array([0.0, mx[1], 0.0]), np.array([0.0, 1.0, 0.0])),
        (np.array([0.0, 0.0, mn[2]]), np.array([0.0, 0.0, -1.0])),
        (np.array([0.0, 0.0, mx[2]]), np.array([0.0, 0.0, 1.0])),
    ]
    return build_polyhedron_from_planes(planes)


def distance_to_box(x: ArrayLike, min_corner: ArrayLike, max_corner: ArrayLike) -> float:
    """Euclidean distance from point x to a closed axis-aligned box."""
    x_arr = _as_vec3("x", x)
    mn = _as_vec3("min_corner", min_corner)
    mx = _as_vec3("max_corner", max_corner)

    delta = np.maximum(np.maximum(mn - x_arr, 0.0), x_arr - mx)
    return float(np.linalg.norm(delta))


def closest_point_on_box_boundary(
    x: ArrayLike,
    min_corner: ArrayLike,
    max_corner: ArrayLike,
) -> FloatArray:
    """Closest point on boundary of an axis-aligned box."""
    x_arr = _as_vec3("x", x)
    mn = _as_vec3("min_corner", min_corner)
    mx = _as_vec3("max_corner", max_corner)

    c = np.clip(x_arr, mn, mx)
    outside = np.any(x_arr < mn) or np.any(x_arr > mx)
    if outside:
        return c

    on_boundary = np.any(np.isclose(c, mn)) or np.any(np.isclose(c, mx))
    if on_boundary:
        return c

    dist_to_min = c - mn
    dist_to_max = mx - c
    stacked = np.stack([dist_to_min, dist_to_max], axis=0)
    side, axis = np.unravel_index(np.argmin(stacked), stacked.shape)

    y = c.copy()
    if side == 0:
        y[axis] = mn[axis]
    else:
        y[axis] = mx[axis]
    return y
