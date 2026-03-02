from __future__ import annotations

import numpy as np

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.visualization import build_polyhedron_mesh, compute_polyhedron_vertices


def test_mesh_reconstruction_for_axis_aligned_box() -> None:
    poly = orient_normals(
        make_axis_aligned_box(
            min_corner=np.array([-1.0, -1.0, -1.0], dtype=float),
            max_corner=np.array([1.0, 1.0, 1.0], dtype=float),
        ),
        interior_point=np.array([0.0, 0.0, 0.0], dtype=float),
    )

    verts = compute_polyhedron_vertices(poly)
    assert verts.shape == (8, 3)

    # Vertices of unit box are exactly +-1 along each axis.
    assert np.all(np.isin(np.abs(verts), [1.0]))

    _, faces = build_polyhedron_mesh(poly)
    assert len(faces) == 6
    assert all(face.shape[0] == 4 for face in faces)

    # Every reconstructed face vertex must satisfy all half-space inequalities.
    for face in faces:
        for p in face:
            d = poly.signed_distances(p)
            assert np.all(d <= 1e-7)
