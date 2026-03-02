import argparse
import itertools

import numpy as np

from wop.geometry import build_polyhedron_from_planes, orient_normals
from wop.visualization import build_polyhedron_mesh, plot_polyhedron


def _base_normals() -> list[np.ndarray]:
    normals: list[np.ndarray] = []

    # 6 axis directions.
    for axis in range(3):
        for sign in (-1.0, 1.0):
            n = np.zeros(3, dtype=float)
            n[axis] = sign
            normals.append(n)

    # 12 edge-diagonal directions.
    edge_set: set[tuple[float, float, float]] = set()
    for sx, sy in itertools.product((-1.0, 1.0), repeat=2):
        base = np.array([sx, sy, 0.0], dtype=float)
        for perm in set(itertools.permutations(base.tolist(), 3)):
            edge_set.add(tuple(float(v) for v in perm))
    normals.extend(np.array(v, dtype=float) for v in sorted(edge_set))

    # 8 body-diagonal directions.
    for sx, sy, sz in itertools.product((-1.0, 1.0), repeat=3):
        normals.append(np.array([sx, sy, sz], dtype=float))

    return normals


def build_complex_polyhedron_planes(
    b_axis: float = 0.9,
    b_edge_diag: float = 1.0,
    b_body_diag: float = 1.1,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Symmetric 26-plane convex polyhedron (all plane families active)."""
    planes: list[tuple[np.ndarray, np.ndarray]] = []
    for n in _base_normals():
        nu = n / np.linalg.norm(n)
        nonzero_count = int(np.count_nonzero(np.abs(n) > 1e-12))
        if nonzero_count == 1:
            b = b_axis
        elif nonzero_count == 2:
            b = b_edge_diag
        else:
            b = b_body_diag
        planes.append((b * nu, nu))
    return planes


def build_asymmetric_polyhedron_planes(
    base: float = 1.0,
    anisotropy: float = 0.22,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Asymmetric convex polyhedron built from same 26 normals.

    Offsets are varied by direction, but kept positive so the origin remains
    in the interior and the shape stays bounded for visualization checks.
    """
    planes: list[tuple[np.ndarray, np.ndarray]] = []
    for n in _base_normals():
        nu = n / np.linalg.norm(n)
        skew = (
            0.80 * nu[0]
            - 0.35 * nu[1]
            + 0.55 * nu[2]
            + 0.40 * nu[0] * nu[1]
            - 0.25 * nu[1] * nu[2]
            + 0.30 * nu[0] * nu[2]
        )
        b = base + anisotropy * skew
        b = max(0.55, float(b))
        planes.append((b * nu, nu))
    return planes


def build_unit_cube_planes() -> list[tuple[np.ndarray, np.ndarray]]:
    """Unit cube [0,1]^3 represented by six boundary planes."""
    return [
        (np.array([0.0, 0.0, 0.0], dtype=float), np.array([-1.0, 0.0, 0.0], dtype=float)),
        (np.array([1.0, 0.0, 0.0], dtype=float), np.array([1.0, 0.0, 0.0], dtype=float)),
        (np.array([0.0, 0.0, 0.0], dtype=float), np.array([0.0, -1.0, 0.0], dtype=float)),
        (np.array([0.0, 1.0, 0.0], dtype=float), np.array([0.0, 1.0, 0.0], dtype=float)),
        (np.array([0.0, 0.0, 0.0], dtype=float), np.array([0.0, 0.0, -1.0], dtype=float)),
        (np.array([0.0, 0.0, 1.0], dtype=float), np.array([0.0, 0.0, 1.0], dtype=float)),
    ]


def _run_case(example: str) -> None:
    if example == "symmetric":
        planes = build_complex_polyhedron_planes()
        title = "Симметричный многогранник"
        face_color = "#74B9FF"
        edge_color = "#0B3D91"
    elif example == "asymmetric":
        planes = build_asymmetric_polyhedron_planes()
        title = "Ассиметричный многогранник"
        face_color = "#8DE1A6"
        edge_color = "#0E5A39"
    elif example == "cube":
        planes = build_unit_cube_planes()
        title = "[0,1]^3"
        face_color = "#FFD166"
        edge_color = "#7A4E00"
    else:
        raise ValueError(f"Unsupported example: {example}")

    poly = orient_normals(build_polyhedron_from_planes(planes), interior_point=np.zeros(3))
    vertices, faces = build_polyhedron_mesh(poly)
    print(
        f"example={example} | "
        f"faces(input planes)={poly.num_faces}, "
        f"vertices={vertices.shape[0]}, visible faces={len(faces)}"
    )

    plot_polyhedron(
        poly,
        title=title,
        face_color=face_color,
        edge_color=edge_color,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Visual test scenes for convex polyhedra.")
    parser.add_argument(
        "--example",
        choices=("symmetric", "asymmetric", "cube"),
        default="asymmetric",
        help="Which polyhedron to render.",
    )
    args = parser.parse_args()
    _run_case(args.example)


if __name__ == "__main__":
    main()
