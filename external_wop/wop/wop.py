from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .estimation import EstimateResult, TrajectoryResult, estimate_from_trajectories
from .geometry import Polyhedron
from .sampling import sample_tangent_direction

FloatArray = NDArray[np.float64]
BoundaryFunc = Callable[[FloatArray, Optional[int]], float]


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


def _sample_plane_radius(
    d: float,
    rng: np.random.Generator,
    min_one_minus_alpha: float,
) -> float:
    """Sample Poisson-kernel radial displacement on a plane at distance d."""
    if d <= 0.0:
        raise ValueError("d must be positive.")
    if min_one_minus_alpha <= 0.0:
        raise ValueError("min_one_minus_alpha must be positive.")

    alpha = float(rng.random())
    one_minus_alpha = max(1.0 - alpha, min_one_minus_alpha)
    rho_sq = d * d * (1.0 / (one_minus_alpha * one_minus_alpha) - 1.0)
    return float(np.sqrt(max(rho_sq, 0.0)))


def trace_wop_trajectory(
    poly: Polyhedron,
    x0: ArrayLike,
    boundary_f: BoundaryFunc,
    rng: np.random.Generator,
    eps_in: float,
    eps_plane: float,
    min_abs_denom: float = 1e-14,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
    r_max: float | None = None,
) -> TrajectoryResult:

    if max_steps <= 0:
        raise ValueError("max_steps must be positive.")
    if eps_in < 0.0 or eps_plane < 0.0:
        raise ValueError("eps_in and eps_plane must be non-negative.")
    if min_abs_denom <= 0.0:
        raise ValueError("min_abs_denom must be positive.")

    x = _as_vec3("x0", x0).copy()

    try:
        i = poly.closest_outside_face_index(x, eps_in=eps_in)
    except ValueError:
        if poly.is_inside_or_on(x, eps_in=eps_in):
            d = poly.signed_distances(x)
            i0 = int(np.argmin(np.abs(d)))
            return TrajectoryResult(float(boundary_f(x, i0)), steps=0, status="hit_face")
        raise ValueError("x0 must belong to the exterior domain.")

    for step in range(1, max_steps + 1):
        if not np.all(np.isfinite(x)):
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")

        if r_max is not None and float(np.linalg.norm(x)) >= r_max:
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")

        nu_i = poly.nu[i]
        b_i = float(poly.b[i])
        d_i = float(np.dot(nu_i, x) - b_i)

        if d_i <= eps_in:
            d_now = poly.signed_distances(x)
            outside_mask_now = d_now > eps_in
            if not np.any(outside_mask_now):
                i_hit = int(np.argmin(np.abs(d_now)))
                return TrajectoryResult(
                    value=float(boundary_f(x, i_hit)),
                    steps=step - 1,
                    status="hit_face",
                )
            i = int(np.argmin(np.where(outside_mask_now, d_now, np.inf)))
            nu_i = poly.nu[i]
            b_i = float(poly.b[i])
            d_i = float(np.dot(nu_i, x) - b_i)

        try:
            w = sample_tangent_direction(
                nu=nu_i,
                rng=rng,
                min_norm=min_abs_denom,
            )
        except RuntimeError:
            return TrajectoryResult(value=float(u_inf), steps=step, status="timeout")

        p = x - d_i * nu_i
        rho = _sample_plane_radius(
            d=d_i,
            rng=rng,
            min_one_minus_alpha=min_abs_denom,
        )
        y = p + rho * w

        # Project back to plane i to reduce floating-point drift.
        y = y - float(np.dot(nu_i, y) - b_i) * nu_i
        if not np.all(np.isfinite(y)):
            return TrajectoryResult(value=float(u_inf), steps=step, status="escaped")

        d = poly.signed_distances(y)
        outside_mask = d > eps_in

        hit = abs(float(d[i])) <= eps_plane and bool(np.all(d <= eps_in))
        if hit:
            return TrajectoryResult(value=float(boundary_f(y, i)), steps=step, status="hit_face")

        if not np.any(outside_mask):
            i_hit = int(np.argmin(np.abs(d)))
            return TrajectoryResult(value=float(boundary_f(y, i_hit)), steps=step, status="hit_face")

        i = int(np.argmin(np.where(outside_mask, d, np.inf)))
        x = y

    return TrajectoryResult(value=float(u_inf), steps=max_steps, status="timeout")


def estimate_wop(
    poly: Polyhedron,
    x0: ArrayLike,
    boundary_f: BoundaryFunc,
    n_paths: int,
    rng: np.random.Generator,
    eps_in: float | None = None,
    eps_plane: float | None = None,
    min_abs_denom: float = 1e-14,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
    r_max: float | None = None,
) -> EstimateResult:

    L = poly.characteristic_length()
    eps_in_eff = 1e-12 * L if eps_in is None else float(eps_in)
    eps_plane_eff = 1e-12 * L if eps_plane is None else float(eps_plane)

    def trace_once() -> TrajectoryResult:
        return trace_wop_trajectory(
            poly=poly,
            x0=x0,
            boundary_f=boundary_f,
            rng=rng,
            eps_in=eps_in_eff,
            eps_plane=eps_plane_eff,
            min_abs_denom=min_abs_denom,
            max_steps=max_steps,
            u_inf=u_inf,
            r_max=r_max,
        )

    return estimate_from_trajectories(n_paths=n_paths, trace_once=trace_once)
