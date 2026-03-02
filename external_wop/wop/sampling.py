from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


def sample_unit_sphere(rng: np.random.Generator) -> FloatArray:
    """Sample isotropic direction on S^2."""
    while True:
        w = rng.normal(size=3)
        n = float(np.linalg.norm(w))
        if n > 0.0:
            return (w / n).astype(float)


def sample_tangent_direction(
    nu: ArrayLike,
    rng: np.random.Generator,
    min_norm: float = 1e-14,
    max_resample: int = 10000,
) -> FloatArray:
    """Sample isotropic unit direction in a plane orthogonal to nu.

    Sampling rule:
      1) draw g ~ N(0, I3)
      2) project to tangent plane: w_raw = g - (g·nu) nu
      3) normalize w_raw
    """
    if min_norm <= 0.0:
        raise ValueError("min_norm must be positive.")
    if max_resample <= 0:
        raise ValueError("max_resample must be positive.")

    nu_arr = _as_vec3("nu", nu)
    nu_norm = float(np.linalg.norm(nu_arr))
    if nu_norm == 0.0:
        raise ValueError("nu must have non-zero norm.")
    nu_unit = (nu_arr / nu_norm).astype(float)

    for _ in range(max_resample):
        g = rng.normal(size=3)
        w_raw = g - float(np.dot(g, nu_unit)) * nu_unit
        n = float(np.linalg.norm(w_raw))
        if not np.isfinite(n) or n <= min_norm:
            continue
        return (w_raw / n).astype(float)

    raise RuntimeError("Failed to sample stable tangent direction.")


def sample_hit_on_plane_from_point(
    x: ArrayLike,
    nu: ArrayLike,
    b: float,
    rng: np.random.Generator,
    min_abs_denom: float = 1e-14,
    max_resample: int = 10000,
) -> tuple[FloatArray, FloatArray, float]:
    """Sample first hit point of an oriented random ray with a plane.

    The plane is ``dot(nu, y)=b`` with unit normal ``nu``.
    Input ``x`` must satisfy ``d=dot(nu,x)-b > 0`` (point on outer side).

    Returns
    -------
    (y, omega, t)
        Hit point, sampled direction, and ray parameter where
        ``y = x + t*omega``.
    """
    x_arr = _as_vec3("x", x)
    nu_arr = _as_vec3("nu", nu)

    d = float(np.dot(nu_arr, x_arr) - b)
    if d <= 0.0:
        raise ValueError("x must satisfy dot(nu,x)-b > 0 for plane sampling.")

    for _ in range(max_resample):
        omega = sample_unit_sphere(rng)
        if float(np.dot(nu_arr, omega)) > 0.0:
            omega = -omega

        denom = float(np.dot(nu_arr, omega))
        if abs(denom) < min_abs_denom:
            continue

        t = -d / denom
        if t <= 0.0 or not np.isfinite(t):
            continue

        y = x_arr + t * omega
        if not np.all(np.isfinite(y)):
            continue
        return y.astype(float), omega.astype(float), float(t)

    raise RuntimeError("Failed to sample stable ray-plane intersection.")
