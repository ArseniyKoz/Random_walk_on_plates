from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Callable, Literal, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .estimation import EstimateResult, TrajectoryResult, estimate_from_trajectories
from .geometry import Polyhedron
from .sampling import sample_tangent_direction

FloatArray = NDArray[np.float64]
BoundaryFunc = Callable[[FloatArray, Optional[int]], float]
RMaxMode = Literal["escape", "project"]


@dataclass(frozen=True)
class _DistanceScanInfo:
    any_outside: bool
    all_inside: bool
    argmin_outside: int | None
    argmin_abs: int
    hit_target_face: bool


@dataclass(frozen=True)
class _RMaxProjection:
    enabled: bool
    center: FloatArray
    radius: float
    radius_sq: float


def _as_vec3(name: str, x: ArrayLike) -> FloatArray:
    arr = np.asarray(x, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}.")
    return arr


def _scan_signed_distances(
    d: FloatArray,
    eps_in: float,
    target_idx: int,
    eps_plane: float,
) -> _DistanceScanInfo:
    if d.ndim != 1 or d.size == 0:
        raise ValueError("d must have shape (M,) with M>0.")
    if target_idx < 0 or target_idx >= d.size:
        raise ValueError("target_idx is out of range.")

    any_outside = False
    all_inside = True
    argmin_outside: int | None = None
    best_outside = float(np.inf)

    argmin_abs = 0
    best_abs = abs(float(d[0]))

    for idx in range(int(d.size)):
        value = float(d[idx])
        if value > eps_in:
            any_outside = True
            all_inside = False
            if value < best_outside:
                best_outside = value
                argmin_outside = idx
        elif not np.isfinite(value):
            all_inside = False

        abs_value = abs(value)
        if abs_value < best_abs:
            best_abs = abs_value
            argmin_abs = idx

    hit_target_face = abs(float(d[target_idx])) <= eps_plane and all_inside
    return _DistanceScanInfo(
        any_outside=any_outside,
        all_inside=all_inside,
        argmin_outside=argmin_outside,
        argmin_abs=argmin_abs,
        hit_target_face=hit_target_face,
    )


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
    linear_tol: float,
) -> FloatArray | None:
    det_a = _det3_rows(nu0, nu1, nu2)
    if abs(det_a) <= linear_tol:
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


def _compute_polyhedron_vertices(
    poly: Polyhedron,
    linear_tol: float = 1e-12,
    inside_tol: float = 1e-9,
    dedup_tol: float = 1e-8,
) -> FloatArray:
    if linear_tol <= 0.0 or inside_tol <= 0.0 or dedup_tol <= 0.0:
        raise ValueError("All tolerances must be positive.")

    m = poly.num_faces
    if m < 4:
        return np.empty((0, 3), dtype=float)

    vertices: list[FloatArray] = []
    dedup_tol_sq = dedup_tol * dedup_tol

    for i, j, k in itertools.combinations(range(m), 3):
        x = _intersect_three_planes(
            nu0=poly.nu[i],
            b0=float(poly.b[i]),
            nu1=poly.nu[j],
            b1=float(poly.b[j]),
            nu2=poly.nu[k],
            b2=float(poly.b[k]),
            linear_tol=linear_tol,
        )
        if x is None:
            continue

        d = poly.signed_distances(x)
        if not bool(np.all(d <= inside_tol)):
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


def _resolve_r_max_projection(
    poly: Polyhedron,
    x0: ArrayLike,
    r_max: float | None,
    r_max_mode: RMaxMode,
    r_max_factor: float,
) -> _RMaxProjection:
    if r_max_mode not in ("escape", "project"):
        raise ValueError("r_max_mode must be one of {'escape', 'project'}.")
    if r_max is not None and float(r_max) <= 0.0:
        raise ValueError("r_max must be positive when provided.")

    if r_max_mode == "escape":
        return _RMaxProjection(
            enabled=False,
            center=np.zeros(3, dtype=float),
            radius=0.0,
            radius_sq=0.0,
        )

    if r_max_factor <= 1.0 or not np.isfinite(r_max_factor):
        raise ValueError("r_max_factor must be finite and greater than 1.0.")

    center, poly_radius = _compute_polyhedron_bounding_sphere(poly)
    radius = float(r_max) if r_max is not None else float(r_max_factor * poly_radius)
    if r_max is None:
        x0_arr = _as_vec3("x0", x0)
        dist_x0 = float(np.linalg.norm(x0_arr - center))
        eps_floor = 1e-12 * max(1.0, poly.characteristic_length())
        radius = max(radius, 1.01 * dist_x0, eps_floor)

    if radius <= 0.0 or not np.isfinite(radius):
        raise ValueError("Effective r_max must be positive and finite.")

    return _RMaxProjection(
        enabled=True,
        center=center,
        radius=radius,
        radius_sq=radius * radius,
    )


def _project_to_sphere(
    x: FloatArray,
    center: FloatArray,
    radius: float,
    min_norm: float,
) -> FloatArray:
    if radius <= 0.0 or not np.isfinite(radius):
        raise ValueError("radius must be positive and finite.")
    if min_norm < 0.0 or not np.isfinite(min_norm):
        raise ValueError("min_norm must be finite and non-negative.")

    dx = x - center
    norm_dx = float(np.linalg.norm(dx))
    if not np.isfinite(norm_dx) or norm_dx <= min_norm:
        return center + np.array([radius, 0.0, 0.0], dtype=float)
    return center + (radius / norm_dx) * dx


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
    r_max_mode: RMaxMode = "escape",
    r_max_factor: float = 3.0,
) -> TrajectoryResult:

    if max_steps <= 0:
        raise ValueError("max_steps must be positive.")
    if eps_in < 0.0 or eps_plane < 0.0:
        raise ValueError("eps_in and eps_plane must be non-negative.")
    if min_abs_denom <= 0.0:
        raise ValueError("min_abs_denom must be positive.")

    x = _as_vec3("x0", x0).copy()
    nu = poly.nu
    b = poly.b
    r_max_sq_escape = None
    if r_max_mode == "escape" and r_max is not None:
        r_max_sq_escape = float(r_max) * float(r_max)
    r_max_cfg = _resolve_r_max_projection(
        poly=poly,
        x0=x,
        r_max=r_max,
        r_max_mode=r_max_mode,
        r_max_factor=r_max_factor,
    )

    try:
        i = poly.closest_outside_face_index(x, eps_in=eps_in)
    except ValueError:
        if poly.is_inside_or_on(x, eps_in=eps_in):
            d = poly.signed_distances(x)
            i0 = _scan_signed_distances(d=d, eps_in=eps_in, target_idx=0, eps_plane=eps_plane).argmin_abs
            return TrajectoryResult(float(boundary_f(x, i0)), steps=0, status="hit_face")
        raise ValueError("x0 must belong to the exterior domain.")

    for step in range(1, max_steps + 1):
        if not np.all(np.isfinite(x)):
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")

        if r_max_cfg.enabled and float(np.dot(x - r_max_cfg.center, x - r_max_cfg.center)) >= r_max_cfg.radius_sq:
            x = _project_to_sphere(
                x=x,
                center=r_max_cfg.center,
                radius=r_max_cfg.radius,
                min_norm=min_abs_denom,
            )
            if not np.all(np.isfinite(x)):
                return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")
            try:
                i = poly.closest_outside_face_index(x, eps_in=eps_in)
            except ValueError:
                if poly.is_inside_or_on(x, eps_in=eps_in):
                    d_proj = poly.signed_distances(x)
                    i_hit = _scan_signed_distances(
                        d=d_proj,
                        eps_in=eps_in,
                        target_idx=0,
                        eps_plane=eps_plane,
                    ).argmin_abs
                    return TrajectoryResult(
                        value=float(boundary_f(x, i_hit)),
                        steps=step - 1,
                        status="hit_face",
                    )
                raise ValueError("Projected point must belong to the exterior domain.")
        elif r_max_sq_escape is not None and float(np.dot(x, x)) >= r_max_sq_escape:
            return TrajectoryResult(value=float(u_inf), steps=step - 1, status="escaped")

        nu_i = nu[i]
        b_i = float(b[i])
        d_i = float(np.dot(nu_i, x) - b_i)

        if d_i <= eps_in:
            d_now = poly.signed_distances(x)
            scan_now = _scan_signed_distances(d=d_now, eps_in=eps_in, target_idx=i, eps_plane=eps_plane)
            if not scan_now.any_outside:
                i_hit = scan_now.argmin_abs
                return TrajectoryResult(
                    value=float(boundary_f(x, i_hit)),
                    steps=step - 1,
                    status="hit_face",
                )
            i = int(scan_now.argmin_outside)  # scan_now.any_outside guarantees non-None.
            nu_i = nu[i]
            b_i = float(b[i])
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
        scan = _scan_signed_distances(d=d, eps_in=eps_in, target_idx=i, eps_plane=eps_plane)
        if scan.hit_target_face:
            return TrajectoryResult(value=float(boundary_f(y, i)), steps=step, status="hit_face")

        if not scan.any_outside:
            i_hit = scan.argmin_abs
            return TrajectoryResult(value=float(boundary_f(y, i_hit)), steps=step, status="hit_face")

        i = int(scan.argmin_outside)  # scan.any_outside guarantees non-None.
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
    r_max_mode: RMaxMode = "escape",
    r_max_factor: float = 3.0,
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
            r_max_mode=r_max_mode,
            r_max_factor=r_max_factor,
        )

    return estimate_from_trajectories(n_paths=n_paths, trace_once=trace_once)
