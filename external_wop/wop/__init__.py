from .geometry import (
    Polyhedron,
    build_polyhedron_from_planes,
    closest_point_on_box_boundary,
    distance_to_box,
    make_axis_aligned_box,
    orient_normals,
)
from .wop import EstimateResult, TrajectoryResult, estimate_wop, trace_wop_trajectory
from .wos_box import estimate_wos_box, trace_wos_box_trajectory
from .visualization import (
    build_face_polygons,
    build_polyhedron_mesh,
    compute_polyhedron_vertices,
    plot_polyhedron,
)

__all__ = [
    "Polyhedron",
    "build_polyhedron_from_planes",
    "closest_point_on_box_boundary",
    "distance_to_box",
    "make_axis_aligned_box",
    "orient_normals",
    "TrajectoryResult",
    "EstimateResult",
    "trace_wop_trajectory",
    "estimate_wop",
    "trace_wos_box_trajectory",
    "estimate_wos_box",
    "compute_polyhedron_vertices",
    "build_face_polygons",
    "build_polyhedron_mesh",
    "plot_polyhedron",
]
