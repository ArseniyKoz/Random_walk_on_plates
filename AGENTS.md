Always first write tests then implement features.
Main language of this repository is C++.

# Repository Guidelines

## Project Structure & Module Organization
Primary code lives in `external_wop_cpp/` (C++20 implementation):
- `src/`: core geometry, sampling, solver, estimation, RNG.
- `include/wop/`: public headers.
- `app/main.cpp`: command-line entry point (`wop_cli`).
- `tests/`: C++ smoke/stats/pipeline tests.

Python package in `external_wop/` is secondary/reference:
- `geometry.py`, `sampling.py`, `wop.py`, `wos_box.py`: core math and Monte Carlo algorithms.
- `estimation.py`: trajectory aggregation/statistics.
- `cli.py`: command-line entry point.

Supporting materials:
- `external_wop/docs/`: derivations and usage notes.
- `external_wop/try/`: exploratory scripts (not part of the package API).
- `external_wop/pyproject.toml`: dependencies, packaging, and pytest config.

## Build, Test, and Development Commands
Run C++ commands from `external_wop_cpp/`.

```bash
/home/admin/Random_walk_on_plates/.venv_test/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWOP_BUILD_TESTS=ON
/home/admin/Random_walk_on_plates/.venv_test/bin/cmake --build build -j
```
Configure and build C++ CLI + tests.

```bash
./build/wop_cli --example box --x0 "3 0 0" --n 50000 --json
```
Run C++ Monte Carlo solver.

```bash
/home/admin/Random_walk_on_plates/.venv_test/bin/ctest --test-dir build --output-on-failure
```
Run C++ tests.

Python (secondary) commands from `external_wop/`.

```bash
python -m pip install -e .[test]
```
Install package in editable mode with test dependencies.

```bash
python -m pip install -e .[viz]
```
Install optional plotting dependencies.

```bash
pytest -q
```
Run the full test suite (configured via `tool.pytest.ini_options`).

```bash
python -m wop.cli --example box --x0 "3 0 0" --n 50000
```
Run a reproducible CLI example and inspect estimator diagnostics.

## Coding Style & Naming Conventions
Target C++20 for production code (`external_wop_cpp/`): prefer clear ownership, const-correctness, and predictable numeric behavior.

Target Python `>=3.10` for reference tooling. Follow PEP 8 with 4-space indentation and explicit type hints (including NumPy types where practical). Use:
- `snake_case` for functions/variables/modules
- `PascalCase` for classes/dataclasses
- clear, validated numeric inputs (shape/type checks for vectors and parameters)

Keep primary algorithmic logic in C++ modules (`external_wop_cpp/src`). Keep CLI/parsing concerns in `app/main.cpp`.

## Testing Guidelines
Use C++ tests in `external_wop_cpp/tests` (via `ctest`) as the primary quality gate. Favor deterministic tests with fixed seeds.

Python tests in `external_wop/tests/` remain useful for parity/regression checks.

## Commit & Pull Request Guidelines
Git history is not included in this workspace snapshot, so no project-specific commit pattern can be inferred here. Use Conventional Commit style by default (for example, `feat: add r_max guard in trajectory tracing`, `fix: tighten plane-hit tolerance check`).

PRs should include:
- concise problem/solution summary
- linked issue (if available)
- test evidence (`pytest -q` output)
- notes on numerical impact and, for visualization changes, a screenshot or plot.
