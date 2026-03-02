from __future__ import annotations

import itertools
from typing import Any, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .geometry import Polyhedron

FloatArray = NDArray[np.float64]


def _normalize(v: FloatArray) -> FloatArray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("Cannot normalize zero vector.")
    return (v / n).astype(float)


def _deduplicate_points(points: FloatArray, tol: float) -> FloatArray:
    if points.size == 0:
        return np.empty((0, 3), dtype=float)

    unique: list[FloatArray] = []
    for p in points:
        if not unique:
            unique.append(np.asarray(p, dtype=float))
            continue

        min_dist = min(float(np.linalg.norm(p - q)) for q in unique)
        if min_dist > tol:
            unique.append(np.asarray(p, dtype=float))

    return np.vstack(unique).astype(float)


def _plane_basis(nu: FloatArray) -> tuple[FloatArray, FloatArray]:
    n = _normalize(np.asarray(nu, dtype=float))
    # Pick a reference axis not almost parallel to n.
    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(ref, n))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)

    u = _normalize(np.cross(n, ref))
    v = _normalize(np.cross(n, u))
    return u, v


def _nice_step(value: float) -> float:
    if value <= 0.0 or not np.isfinite(value):
        return 1.0

    exponent = np.floor(np.log10(value))
    fraction = value / (10.0 ** exponent)

    if fraction <= 1.0:
        nice = 1.0
    elif fraction <= 2.0:
        nice = 2.0
    elif fraction <= 5.0:
        nice = 5.0
    else:
        nice = 10.0

    return float(nice * (10.0 ** exponent))


def _auto_tick_values(radius: float, tick_step: Optional[float]) -> FloatArray:
    if tick_step is None:
        tick_step = _nice_step(radius / 3.0)
    if tick_step <= 0.0 or not np.isfinite(tick_step):
        tick_step = 1.0

    max_tick = tick_step * np.floor(radius / tick_step)
    if max_tick <= 0.0:
        return np.array([0.0], dtype=float)

    ticks = np.arange(-max_tick, max_tick + 0.5 * tick_step, tick_step, dtype=float)
    if ticks.size > 21:
        # Keep labels readable in very large scenes.
        stride = int(np.ceil(ticks.size / 21.0))
        ticks = ticks[::stride]
    return ticks.astype(float)


def _draw_stable_axes(
    ax: Any,
    origin: FloatArray,
    radius: float,
    tick_values: Optional[FloatArray],
    tick_step: Optional[float],
    show_tick_points: bool,
) -> None:
    dirs = np.eye(3, dtype=float)
    colors = ("#C1121F", "#2A9D8F", "#1D4ED8")
    labels = ("X", "Y", "Z")

    ticks = _auto_tick_values(radius, tick_step) if tick_values is None else np.asarray(tick_values, dtype=float)
    ticks = ticks[np.isfinite(ticks)]

    axis_extent = 1.05 * radius
    label_offset = 0.07 * radius

    for idx in range(3):
        d = dirs[idx]
        p0 = origin - axis_extent * d
        p1 = origin + axis_extent * d
        ax.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            [p0[2], p1[2]],
            color=colors[idx],
            linewidth=1.6,
        )

        if show_tick_points and ticks.size > 0:
            pts = origin[None, :] + np.outer(ticks, d)
            ax.scatter(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                s=16.0,
                c=colors[idx],
                depthshade=False,
            )

            # Draw compact numeric labels except for zero to avoid clutter.
            offset = label_offset * (dirs[(idx + 1) % 3] + 0.45 * dirs[(idx + 2) % 3])
            for t, p in zip(ticks, pts):
                if abs(float(t)) < 1e-14:
                    continue
                ax.text(
                    float(p[0] + offset[0]),
                    float(p[1] + offset[1]),
                    float(p[2] + offset[2]),
                    f"{t:g}",
                    color=colors[idx],
                    fontsize=8,
                )

        p_label = origin + (axis_extent + label_offset) * d
        ax.text(
            float(p_label[0]),
            float(p_label[1]),
            float(p_label[2]),
            labels[idx],
            color=colors[idx],
            fontsize=10,
            fontweight="bold",
        )


def order_coplanar_points(points: FloatArray, nu: FloatArray) -> FloatArray:
    """Order coplanar points around centroid in a stable angular order."""
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (K,3).")
    if points.shape[0] < 3:
        raise ValueError("At least 3 points are required to form a polygon.")

    c = np.mean(points, axis=0)
    u, v = _plane_basis(nu)

    rel = points - c
    x = rel @ u
    y = rel @ v
    ang = np.arctan2(y, x)
    order = np.argsort(ang)
    return points[order].astype(float)


def compute_polyhedron_vertices(
    poly: Polyhedron,
    linear_tol: float = 1e-12,
    inside_tol: float = 1e-9,
    dedup_tol: float = 1e-8,
) -> FloatArray:
    """Reconstruct polyhedron vertices from plane intersections.

    Vertices are found as intersections of all non-parallel triplets of planes,
    then filtered by half-space constraints ``nu @ x - b <= inside_tol``.
    """
    if linear_tol <= 0.0 or inside_tol <= 0.0 or dedup_tol <= 0.0:
        raise ValueError("All tolerances must be positive.")

    m = poly.num_faces
    if m < 4:
        return np.empty((0, 3), dtype=float)

    vertices: list[FloatArray] = []
    idx = range(m)
    for i, j, k in itertools.combinations(idx, 3):
        a = np.vstack((poly.nu[i], poly.nu[j], poly.nu[k]))
        det = float(np.linalg.det(a))
        if abs(det) <= linear_tol:
            continue

        rhs = np.array([poly.b[i], poly.b[j], poly.b[k]], dtype=float)
        try:
            x = np.linalg.solve(a, rhs)
        except np.linalg.LinAlgError:
            continue

        d = poly.signed_distances(x)
        if np.all(d <= inside_tol):
            vertices.append(np.asarray(x, dtype=float))

    if not vertices:
        return np.empty((0, 3), dtype=float)

    points = np.vstack(vertices).astype(float)
    return _deduplicate_points(points, tol=dedup_tol)


def build_face_polygons(
    poly: Polyhedron,
    vertices: Optional[FloatArray] = None,
    plane_tol: float = 1e-7,
    dedup_tol: float = 1e-8,
) -> list[FloatArray]:
    """Build ordered face polygons from plane representation.

    Returns one polygon per non-empty face. Each polygon is an array ``(K,3)``.
    """
    if plane_tol <= 0.0 or dedup_tol <= 0.0:
        raise ValueError("plane_tol and dedup_tol must be positive.")

    verts = compute_polyhedron_vertices(poly) if vertices is None else np.asarray(vertices, dtype=float)
    if verts.ndim != 2 or (verts.size > 0 and verts.shape[1] != 3):
        raise ValueError("vertices must have shape (N,3).")

    if verts.shape[0] == 0:
        return []

    faces: list[FloatArray] = []
    for i in range(poly.num_faces):
        dist = verts @ poly.nu[i] - poly.b[i]
        mask = np.abs(dist) <= plane_tol
        face_pts = verts[mask]
        if face_pts.shape[0] < 3:
            continue

        face_pts = _deduplicate_points(face_pts, tol=dedup_tol)
        if face_pts.shape[0] < 3:
            continue

        polygon = order_coplanar_points(face_pts, poly.nu[i])
        faces.append(polygon)

    return faces


def build_polyhedron_mesh(
    poly: Polyhedron,
    linear_tol: float = 1e-12,
    inside_tol: float = 1e-9,
    dedup_tol: float = 1e-8,
    plane_tol: float = 1e-7,
) -> tuple[FloatArray, list[FloatArray]]:
    """Return reconstructed ``(vertices, face_polygons)`` mesh."""
    verts = compute_polyhedron_vertices(
        poly=poly,
        linear_tol=linear_tol,
        inside_tol=inside_tol,
        dedup_tol=dedup_tol,
    )
    faces = build_face_polygons(
        poly=poly,
        vertices=verts,
        plane_tol=plane_tol,
        dedup_tol=dedup_tol,
    )
    return verts, faces


def plot_polyhedron(
    poly: Polyhedron,
    ax: Any = None,
    title: str = "Polyhedron",
    face_color: str = "#5FA8D3",
    edge_color: str = "#1D3557",
    alpha: float = 0.35,
    show_vertices: bool = True,
    stable_axes: bool = True,
    axis_origin: Optional[ArrayLike] = None,
    axis_ticks: Optional[ArrayLike] = None,
    axis_tick_step: Optional[float] = None,
    hide_matplotlib_axes: bool = True,
    show_axis_points: bool = True,
    view_elev: Optional[float] = 23.33333333333333,
    view_azim: Optional[float] = -125.0,
    view_roll: Optional[float] = -150.0,
    view_vertical_axis: str = "z",
    show: bool = True,
) -> Any:
    """Plot reconstructed convex polyhedron using matplotlib.

    Matplotlib is imported lazily, so core package dependencies stay minimal.
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:
        raise ImportError(
            "plot_polyhedron requires matplotlib. Install it with: pip install matplotlib"
        ) from exc

    verts, faces = build_polyhedron_mesh(poly)
    if verts.shape[0] == 0 or len(faces) == 0:
        raise ValueError(
            "Could not reconstruct a bounded polyhedron mesh from provided planes."
        )

    created_ax = False
    if ax is None:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
        created_ax = True

    coll = Poly3DCollection(
        faces,
        facecolors=face_color,
        edgecolors=edge_color,
        linewidths=1.0,
        alpha=alpha,
        zsort="average",
    )
    ax.add_collection3d(coll)

    if show_vertices:
        ax.scatter(verts[:, 0], verts[:, 1], verts[:, 2], s=18.0, c=edge_color, depthshade=False)

    mins = np.min(verts, axis=0)
    maxs = np.max(verts, axis=0)
    center = 0.5 * (mins + maxs)
    radius = 0.5 * float(np.max(maxs - mins))
    radius = max(radius, 1.0) * 1.1

    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)

    if view_elev is not None and view_azim is not None:
        # Camera preset:
        # - x projects toward lower-left screen corner
        # - y projects approximately parallel to bottom screen edge.
        try:
            ax.view_init(
                elev=float(view_elev),
                azim=float(view_azim),
                roll=float(view_roll) if view_roll is not None else 0.0,
                vertical_axis=view_vertical_axis,
            )
        except TypeError:
            # Fallback for older Matplotlib signatures.
            ax.view_init(elev=float(view_elev), azim=float(view_azim))

    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1.0, 1.0, 1.0))

    if stable_axes:
        origin = np.zeros(3, dtype=float) if axis_origin is None else np.asarray(axis_origin, dtype=float)
        if origin.shape != (3,):
            raise ValueError("axis_origin must have shape (3,).")
        ticks = None if axis_ticks is None else np.asarray(axis_ticks, dtype=float).ravel()
        _draw_stable_axes(
            ax=ax,
            origin=origin,
            radius=radius,
            tick_values=ticks,
            tick_step=axis_tick_step,
            show_tick_points=show_axis_points,
        )
        if hide_matplotlib_axes:
            # Built-in mplot3d axes reassign visible edges while rotating,
            # which looks like ticks/points jumping across faces.
            ax.set_axis_off()
    else:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

    ax.set_title(title)

    if show and created_ax:
        plt.show()

    return ax
