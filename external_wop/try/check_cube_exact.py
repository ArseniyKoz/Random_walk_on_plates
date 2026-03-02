from __future__ import annotations

import argparse

import numpy as np

from wop.argparse_utils import parse_vec3_arg
from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import estimate_wop


def run_check(
    x0_list: list[np.ndarray],
    n_paths: int,
    seed: int,
    max_steps: int,
    r_max: float | None,
    a: np.ndarray,
    u_num: float | None,
) -> None:
    box_min = np.array([0.0, 0.0, 0.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    interior = 0.5 * (box_min + box_max)

    poly = orient_normals(make_axis_aligned_box(box_min, box_max), interior)
    if not poly.is_inside_or_on(a):
        raise ValueError("Point 'a' must be inside or on the unit cube [0,1]^3.")

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    def boundary_f(y: np.ndarray, _face: int | None) -> float:
        return exact_u(y)

    print("Verification on unit cube [0,1]^3")
    print(f"a (source point inside cube): {a}")
    print(f"n_paths={n_paths}, max_steps={max_steps}, r_max={r_max}, base_seed={seed}")
    print()
    print(
        f"{'x0':>24} | {'J (numeric)':>12} | {'exact':>12} | {'|J-exact|':>12} | "
        f"{'eps':>12} | {'ok':>3} | {'truncated':>10}"
    )
    print("-" * 105)

    for idx, x0 in enumerate(x0_list):
        if poly.is_inside_or_on(x0):
            raise ValueError(f"x0 must be outside cube, got {x0}.")

        rng = np.random.default_rng(seed + idx)
        res = estimate_wop(
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
        err = abs(res.J - exact)
        ok = err <= res.eps
        x0_str = np.array2string(x0, precision=3, floatmode="fixed")

        print(
            f"{x0_str:>24} | {res.J:12.6f} | {exact:12.6f} | {err:12.6f} | "
            f"{res.eps:12.6f} | {('yes' if ok else 'no'):>3} | "
            f"{res.n_truncated:5d}/{res.n_total:<4d}"
        )

        if u_num is not None:
            diff_num_exact = abs(u_num - exact)
            diff_num_j = abs(u_num - res.J)
            print(
                f"{'your_u':>24} | {u_num:12.6f} | {'-':>12} | "
                f"{'|u-exact|=' + format(diff_num_exact, '.6f'):>12} | "
                f"{'|u-J|=' + format(diff_num_j, '.6f'):>12} | {'-':>3} | {'-':>10}"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check WOP numeric result against exact harmonic solution outside unit cube."
    )
    parser.add_argument(
        "--x0",
        type=lambda s: parse_vec3_arg(s, example="2 0.2 0.3"),
        action="append",
        help="Outside start point; can be repeated. Example: --x0 '2 0.2 0.3'",
    )
    parser.add_argument(
        "--a",
        type=lambda s: parse_vec3_arg(s, example="0.2 0.3 0.7"),
        default=parse_vec3_arg("0.2 0.3 0.7", example="0.2 0.3 0.7"),
        help="Point a inside cube for exact solution u(x)=1/||x-a||.",
    )
    parser.add_argument("--n", type=int, default=20000, help="Number of trajectories per x0.")
    parser.add_argument("--seed", type=int, default=12345, help="Base RNG seed.")
    parser.add_argument("--max-steps", type=int, default=200000, help="Max steps per trajectory.")
    parser.add_argument(
        "--r-max",
        type=float,
        default=1e6,
        help="Escape radius (set <=0 to disable).",
    )
    parser.add_argument(
        "--u-num",
        type=float,
        default=None,
        help="Optional: your numeric value to compare with exact and with J.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.n <= 0:
        parser.error("--n must be positive.")
    if args.max_steps <= 0:
        parser.error("--max-steps must be positive.")

    x0_list = args.x0
    if x0_list is None or len(x0_list) == 0:
        x0_list = [
            parse_vec3_arg("2.0 0.2 0.3", example="2 0.2 0.3"),
            parse_vec3_arg("-1.5 0.5 0.5", example="2 0.2 0.3"),
            parse_vec3_arg("0.2 0.4 2.0", example="2 0.2 0.3"),
        ]

    run_check(
        x0_list=x0_list,
        n_paths=args.n,
        seed=args.seed,
        max_steps=args.max_steps,
        r_max=(args.r_max if args.r_max is not None and args.r_max > 0.0 else None),
        a=args.a,
        u_num=args.u_num,
    )


if __name__ == "__main__":
    main()
