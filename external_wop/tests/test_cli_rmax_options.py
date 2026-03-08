from __future__ import annotations

import pytest

from wop.cli import build_parser


def test_cli_parses_rmax_projection_flags() -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "--example",
            "box",
            "--x0",
            "3 0 0",
            "--n",
            "10",
            "--r-max",
            "2.5",
            "--r-max-mode",
            "project",
            "--r-max-factor",
            "4.0",
        ]
    )

    assert args.r_max == pytest.approx(2.5)
    assert args.r_max_mode == "project"
    assert args.r_max_factor == pytest.approx(4.0)


def test_cli_has_expected_rmax_defaults() -> None:
    parser = build_parser()
    args = parser.parse_args([])

    assert args.r_max == pytest.approx(1e6)
    assert args.r_max_mode == "escape"
    assert args.r_max_factor == pytest.approx(3.0)
