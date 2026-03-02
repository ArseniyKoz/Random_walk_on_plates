from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .estimation import EstimateResult, TrajectoryResult, estimate_from_trajectories
from .geometry import closest_point_on_box_boundary, distance_to_box
from .sampling import sample_unit_sphere

FloatArray = NDArray[np.float64]
BoundaryFunc = Callable[[FloatArray, Optional[int]], float]


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


def trace_wos_box_trajectory(
    x0: ArrayLike,
    box_min: ArrayLike,
    box_max: ArrayLike,
    boundary_f: BoundaryFunc,
    rng: np.random.Generator,
    eps: float,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
    r_max: float | None = None,
) -> TrajectoryResult:
    """Trace one walk-on-spheres trajectory for exterior of axis-aligned box."""
    if eps <= 0.0:
        raise ValueError("eps must be positive.")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive.")

    x = _as_vec3("x0", x0).copy()
    mn = _as_vec3("box_min", box_min)
    mx = _as_vec3("box_max", box_max)
    r_max_sq = None if r_max is None else float(r_max) * float(r_max)

    for step in range(1, max_steps + 1):
        if r_max_sq is not None and float(np.dot(x, x)) >= r_max_sq:
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")

        dist = distance_to_box(x, mn, mx)
        if dist <= eps:
            y = closest_point_on_box_boundary(x, mn, mx)
            return TrajectoryResult(value=float(boundary_f(y, None)), steps=step, status="hit_face")

        omega = sample_unit_sphere(rng)
        x = x + dist * omega

    return TrajectoryResult(value=float(u_inf), steps=max_steps, status="timeout")


def estimate_wos_box(
    x0: ArrayLike,
    box_min: ArrayLike,
    box_max: ArrayLike,
    boundary_f: BoundaryFunc,
    n_paths: int,
    rng: np.random.Generator,
    eps: float,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
    r_max: float | None = None,
) -> EstimateResult:
    """Monte Carlo estimator for box exterior via walk-on-spheres."""
    def trace_once() -> TrajectoryResult:
        return trace_wos_box_trajectory(
            x0=x0,
            box_min=box_min,
            box_max=box_max,
            boundary_f=boundary_f,
            rng=rng,
            eps=eps,
            max_steps=max_steps,
            u_inf=u_inf,
            r_max=r_max,
        )

    return estimate_from_trajectories(n_paths=n_paths, trace_once=trace_once)
