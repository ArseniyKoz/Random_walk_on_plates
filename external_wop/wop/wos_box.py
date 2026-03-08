from __future__ import annotations

import itertools
from typing import Callable, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .estimation import EstimateResult, TrajectoryResult, estimate_from_trajectories
from .geometry import Polyhedron, make_axis_aligned_box, orient_normals
from .sampling import sample_unit_sphere

FloatArray = NDArray[np.float64]
BoundaryFunc = Callable[[FloatArray, Optional[int]], float]

_LINEAR_TOL = 1e-12
_INSIDE_TOL = 1e-9
_DEDUP_TOL = 1e-8
_LAMBDA_TOL = 1e-10
_FEAS_TOL = 1e-9
_MIN_ORTH_NORM = 1e-14


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


def _det3_rows(r0: FloatArray, r1: FloatArray, r2: FloatArray) -> float:
    return float(
        r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
        - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
        + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0])
    )


def _intersect_three_planes(
    nu0: FloatArray,
    b0: float,
    nu1: FloatArray,
    b1: float,
    nu2: FloatArray,
    b2: float,
) -> FloatArray | None:
    det_a = _det3_rows(nu0, nu1, nu2)
    if abs(det_a) <= _LINEAR_TOL:
        return None

    det_x = _det3_rows(
        np.array([b0, nu0[1], nu0[2]], dtype=float),
        np.array([b1, nu1[1], nu1[2]], dtype=float),
        np.array([b2, nu2[1], nu2[2]], dtype=float),
    )
    det_y = _det3_rows(
        np.array([nu0[0], b0, nu0[2]], dtype=float),
        np.array([nu1[0], b1, nu1[2]], dtype=float),
        np.array([nu2[0], b2, nu2[2]], dtype=float),
    )
    det_z = _det3_rows(
        np.array([nu0[0], nu0[1], b0], dtype=float),
        np.array([nu1[0], nu1[1], b1], dtype=float),
        np.array([nu2[0], nu2[1], b2], dtype=float),
    )
    x = np.array([det_x / det_a, det_y / det_a, det_z / det_a], dtype=float)
    if not np.all(np.isfinite(x)):
        return None
    return x


def _compute_polyhedron_vertices(poly: Polyhedron) -> FloatArray:
    m = poly.num_faces
    if m < 4:
        return np.empty((0, 3), dtype=float)

    vertices: list[FloatArray] = []
    dedup_tol_sq = _DEDUP_TOL * _DEDUP_TOL
    for i, j, k in itertools.combinations(range(m), 3):
        x = _intersect_three_planes(
            nu0=poly.nu[i],
            b0=float(poly.b[i]),
            nu1=poly.nu[j],
            b1=float(poly.b[j]),
            nu2=poly.nu[k],
            b2=float(poly.b[k]),
        )
        if x is None:
            continue

        d = poly.signed_distances(x)
        if not bool(np.all(d <= _INSIDE_TOL)):
            continue

        duplicate = False
        for v in vertices:
            if float(np.dot(x - v, x - v)) <= dedup_tol_sq:
                duplicate = True
                break
        if not duplicate:
            vertices.append(x)

    if not vertices:
        return np.empty((0, 3), dtype=float)
    return np.vstack(vertices).astype(float)


def _compute_polyhedron_bounding_sphere(poly: Polyhedron) -> tuple[FloatArray, float]:
    vertices = _compute_polyhedron_vertices(poly)
    if vertices.shape[0] == 0:
        return np.zeros(3, dtype=float), max(1.0, poly.characteristic_length())

    v_min = np.min(vertices, axis=0)
    v_max = np.max(vertices, axis=0)
    center = (0.5 * (v_min + v_max)).astype(float)
    radius = float(np.max(np.linalg.norm(vertices - center[None, :], axis=1)))
    eps_floor = 1e-12 * max(1.0, poly.characteristic_length())
    return center, max(radius, eps_floor)


def _project_to_polyhedron_boundary(poly: Polyhedron, x: FloatArray) -> tuple[FloatArray, float]:
    m = poly.num_faces
    best_y: FloatArray | None = None
    best_d2 = float(np.inf)

    def try_active_set(active: tuple[int, ...]) -> None:
        nonlocal best_y, best_d2
        n_active = poly.nu[list(active)]
        rhs = n_active @ x - poly.b[list(active)]
        gram = n_active @ n_active.T

        try:
            lam = np.linalg.solve(gram, rhs)
        except np.linalg.LinAlgError:
            return

        if np.any(lam < -_LAMBDA_TOL):
            return

        y = x - n_active.T @ lam
        if not np.all(np.isfinite(y)):
            return

        d = poly.signed_distances(y)
        if np.any(d > _FEAS_TOL):
            return

        d2 = float(np.dot(x - y, x - y))
        if d2 < best_d2:
            best_d2 = d2
            best_y = y.astype(float)

    for i in range(m):
        try_active_set((i,))
    for i in range(m):
        for j in range(i + 1, m):
            try_active_set((i, j))
    for i in range(m):
        for j in range(i + 1, m):
            for k in range(j + 1, m):
                try_active_set((i, j, k))

    if best_y is None:
        raise RuntimeError("Failed to project point to polyhedron boundary.")
    return best_y, float(np.sqrt(max(best_d2, 0.0)))


def _sample_unit_orthogonal(e: FloatArray, rng: np.random.Generator) -> FloatArray:
    for _ in range(512):
        omega = sample_unit_sphere(rng)
        w_raw = omega - float(np.dot(omega, e)) * e
        n = float(np.linalg.norm(w_raw))
        if n > _MIN_ORTH_NORM:
            return (w_raw / n).astype(float)
    raise RuntimeError("Could not sample orthogonal direction.")


def _sample_far_sphere_step(
    x: FloatArray,
    center: FloatArray,
    rho: float,
    rng: np.random.Generator,
) -> FloatArray:
    dx = x - center
    r = float(np.linalg.norm(dx))
    if not r > rho:
        raise ValueError("Far-sphere step requires |x-center| > rho.")
    e = dx / r

    alpha = float(rng.random())
    denom = r - rho + 2.0 * alpha * rho
    if denom <= 0.0:
        raise RuntimeError("Invalid denominator in far-sphere sampling.")

    q = (r * r - rho * rho) / denom
    z1 = (r * r + rho * rho - q * q) / (2.0 * r * rho)
    z1 = float(np.clip(z1, -1.0, 1.0))

    w = _sample_unit_orthogonal(e, rng)
    sin_part = float(np.sqrt(max(1.0 - z1 * z1, 0.0)))
    direction = z1 * e + sin_part * w
    return (center + rho * direction).astype(float)


def _argmin_abs(values: FloatArray) -> int:
    if values.size == 0:
        raise ValueError("values must be non-empty.")
    return int(np.argmin(np.abs(values)))


def _trace_wos_poly_trajectory(
    poly: Polyhedron,
    x0: ArrayLike,
    boundary_f: BoundaryFunc,
    rng: np.random.Generator,
    delta: float,
    rho_scale: float,
    rho1_scale: float,
    max_steps: int,
    u_inf: float,
) -> TrajectoryResult:
    if delta <= 0.0:
        raise ValueError("delta must be positive.")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive.")
    if rho_scale <= 0.0 or not np.isfinite(rho_scale):
        raise ValueError("rho_scale must be finite and positive.")
    if rho1_scale <= 1.0 or not np.isfinite(rho1_scale):
        raise ValueError("rho1_scale must be finite and greater than 1.0.")

    x = _as_vec3("x0", x0).copy()
    try:
        _ = poly.closest_outside_face_index(x, eps_in=0.0)
    except ValueError:
        if poly.is_inside_or_on(x, eps_in=0.0):
            d = poly.signed_distances(x)
            i_hit = _argmin_abs(d)
            return TrajectoryResult(value=float(boundary_f(x, i_hit)), steps=0, status="hit_face")
        raise ValueError("x0 must belong to the exterior domain.")

    center, base_rho = _compute_polyhedron_bounding_sphere(poly)
    rho = float(rho_scale * base_rho)
    rho1 = float(rho1_scale * rho)
    if not (rho > 0.0 and rho1 > rho):
        raise ValueError("Invalid rho/rho1 configuration.")

    eta = 1.0
    for step in range(1, max_steps + 1):
        if not np.all(np.isfinite(x)) or not np.isfinite(eta):
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="timeout")

        x_star, dist = _project_to_polyhedron_boundary(poly, x)
        if dist <= delta:
            return TrajectoryResult(
                value=float(eta * float(boundary_f(x_star, None))),
                steps=step,
                status="hit_face",
            )

        r = float(np.linalg.norm(x - center))
        if r <= rho1:
            omega = sample_unit_sphere(rng)
            x = x + dist * omega
        else:
            x = _sample_far_sphere_step(x=x, center=center, rho=rho, rng=rng)
            eta *= rho / r

    return TrajectoryResult(value=float(u_inf), steps=max_steps, status="timeout")


def trace_wos_box_trajectory(
    x0: ArrayLike,
    box_min: ArrayLike,
    box_max: ArrayLike,
    boundary_f: BoundaryFunc,
    rng: np.random.Generator,
    delta: float,
    rho_scale: float = 1.0,
    rho1_scale: float = 2.0,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
) -> TrajectoryResult:
    """Lecture-style WoS trajectory for exterior of axis-aligned box."""
    mn = _as_vec3("box_min", box_min)
    mx = _as_vec3("box_max", box_max)
    poly = orient_normals(make_axis_aligned_box(mn, mx), interior_point=0.5 * (mn + mx))
    return _trace_wos_poly_trajectory(
        poly=poly,
        x0=x0,
        boundary_f=boundary_f,
        rng=rng,
        delta=delta,
        rho_scale=rho_scale,
        rho1_scale=rho1_scale,
        max_steps=max_steps,
        u_inf=u_inf,
    )


def estimate_wos_box(
    x0: ArrayLike,
    box_min: ArrayLike,
    box_max: ArrayLike,
    boundary_f: BoundaryFunc,
    n_paths: int,
    rng: np.random.Generator,
    delta: float,
    rho_scale: float = 1.0,
    rho1_scale: float = 2.0,
    max_steps: int = 1_000_000,
    u_inf: float = 0.0,
) -> EstimateResult:
    """Lecture-style Monte Carlo estimator for box exterior."""

    mn = _as_vec3("box_min", box_min)
    mx = _as_vec3("box_max", box_max)
    poly = orient_normals(make_axis_aligned_box(mn, mx), interior_point=0.5 * (mn + mx))

    def trace_once() -> TrajectoryResult:
        return _trace_wos_poly_trajectory(
            poly=poly,
            x0=x0,
            boundary_f=boundary_f,
            rng=rng,
            delta=delta,
            rho_scale=rho_scale,
            rho1_scale=rho1_scale,
            max_steps=max_steps,
            u_inf=u_inf,
        )

    return estimate_from_trajectories(n_paths=n_paths, trace_once=trace_once)
