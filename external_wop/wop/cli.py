from __future__ import annotations

import argparse

import numpy as np

from .argparse_utils import parse_vec3_arg
from .geometry import make_axis_aligned_box, orient_normals
from .wop import estimate_wop


def run_box_example(x0: np.ndarray, n_paths: int, seed: int, max_steps: int, r_max: float | None) -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    interior = np.array([0.0, 0.0, 0.0], dtype=float)

    poly = make_axis_aligned_box(box_min, box_max)
    poly = orient_normals(poly, interior)

    a = np.array([0.2, -0.1, 0.3], dtype=float)

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    def boundary_f(y: np.ndarray, _face: int | None) -> float:
        return exact_u(y)

    rng = np.random.default_rng(seed)
    result = estimate_wop(
        poly=poly,
        x0=x0,
        boundary_f=boundary_f,
        n_paths=n_paths,
        rng=rng,
        max_steps=max_steps,
        u_inf=0.0,
        r_max=r_max,
    )

    exact = exact_u(x0)
    abs_err = abs(result.J - exact)

    print(f"x0: {x0}")
    print(f"n_total: {result.n_total}")
    print(f"n_truncated: {result.n_truncated}")
    print(f"J: {result.J:.10f}")
    print(f"exact: {exact:.10f}")
    print(f"abs_error: {abs_err:.10f}")
    print(f"S2: {result.S2:.10f}")
    print(f"eps: {result.eps:.10f}")
    print(f"mean_steps: {result.mean_steps:.2f}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Walk on Planes estimator in R^3 exterior domain.")
    parser.add_argument("--example", choices=["box"], default="box", help="Built-in example.")
    parser.add_argument(
        "--x0",
        type=lambda s: parse_vec3_arg(s, example="3 0 0"),
        default=parse_vec3_arg("3 0 0", example="3 0 0"),
        help="Start point in exterior domain.",
    )
    parser.add_argument("--n", type=int, default=50000, help="Number of Monte Carlo trajectories.")
    parser.add_argument("--seed", type=int, default=12345, help="RNG seed.")
    parser.add_argument("--max-steps", type=int, default=1_000_000, help="Maximum steps per trajectory.")
    parser.add_argument(
        "--r-max",
        type=float,
        default=1e6,
        help="Escape radius used to approximate infinity boundary (set to 0 or negative to disable).",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.n <= 0:
        parser.error("--n must be positive")

    if args.example == "box":
        run_box_example(
            x0=args.x0,
            n_paths=args.n,
            seed=args.seed,
            max_steps=args.max_steps,
            r_max=(args.r_max if args.r_max is not None and args.r_max > 0.0 else None),
        )
    else:
        parser.error(f"Unsupported example: {args.example}")


if __name__ == "__main__":
    main()
