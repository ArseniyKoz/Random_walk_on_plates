# Walk on Planes (WOP) for Exterior Dirichlet Problem in R^3

This project implements a Monte Carlo solver for the exterior Dirichlet problem for a convex polyhedron with planar faces, using the Walk on Planes (WOP) algorithm.

## Features

- Polyhedron geometry stored only as face planes `(p_i, nu_i)` or internal `(nu, b)`.
- Reproducible random sampling via `numpy.random.Generator`.
- WOP trajectory simulation with face-hit stopping criterion.
- WOP plane step with tangent-direction sampling via projected Gaussian vector.
- Reference Walk on Spheres (WOS) implementation for axis-aligned box.
- Polyhedron mesh reconstruction and 3D plotting from plane data only.
- Tests for plane sampler, exact harmonic solution, and WOP-vs-WOS cross-check.

## Install

```bash
pip install -e .[test]
```

Optional visualization dependency:

```bash
pip install -e .[viz]
```

## Run Tests

```bash
pytest -q
```

## C++ Build (Core + CLI)

The C++20 port is now separated from Python package and located in sibling directory `../external_wop_cpp/`.

Configure and build:

```bash
cd ../external_wop_cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Run the C++ CLI:

```bash
./build/wop_cli --example box --x0 "3 0 0" --n 50000
```

Run the C++ CLI from YAML config:

```bash
./build/Release/wop_cli.exe --config examples/box_wop.yaml --json
./build/Release/wop_cli.exe --config examples/box_wos.yaml --json
```

Ready-to-run YAML examples live in `../external_wop_cpp/examples/`.
Each config includes a top-level `command:` field with the intended launch command.
Config files now select built-in boundary functions from C++ code instead of
parsing free-form formulas. The specialized `--example box` path is still
faster than the general config path.

Machine-readable output for regression checks:

```bash
./build/wop_cli --example box --x0 "3 0 0" --n 50000 --json
```

Run C++ tests:

```bash
ctest --test-dir build --output-on-failure
```

Optional Python parity test against C++ CLI:

```bash
WOP_CPP_CLI=../external_wop_cpp/build/wop_cli pytest -q tests/test_cpp_cli_parity.py
```

## CLI Example

```bash
python -m wop.cli --example box --x0 "3 0 0" --n 50000
```

Projection mode with auto-`r_max` from polyhedron geometry:

```bash
python -m wop.cli --example box --x0 "3 0 0" --n 50000 --r-max 0 --r-max-mode project --r-max-factor 3.0
```

Expected output includes Monte Carlo `J`, `S2`, `eps = 3*sqrt(S2/N)`, exact value for built-in example, and trajectory diagnostics including total and truncated walks.

## Visual Check of Polyhedron

```python
import numpy as np
from wop.geometry import build_polyhedron_from_planes, orient_normals
from wop.visualization import plot_polyhedron

planes = [
    (np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])),
    (np.array([-1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])),
    (np.array([0.0, 1.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    (np.array([0.0, -1.0, 0.0]), np.array([0.0, -1.0, 0.0])),
    (np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 1.0])),
    (np.array([0.0, 0.0, -1.0]), np.array([0.0, 0.0, -1.0])),
]
poly = orient_normals(build_polyhedron_from_planes(planes), interior_point=np.zeros(3))
plot_polyhedron(poly, title="Input Polyhedron")
```

By default, `plot_polyhedron` now draws stable custom coordinate axes and hides
Matplotlib's dynamic 3D axes (those often \"jump\" between cube edges while rotating).
To restore old behavior, pass `stable_axes=False`.

## Project Layout

- `wop/geometry.py`: polyhedron representation and geometry utilities.
- `wop/sampling.py`: isotropic direction and ray-plane hit sampling.
- `wop/wop.py`: WOP trajectory and Monte Carlo estimator.
- `wop/wos_box.py`: WOS reference estimator for axis-aligned box.
- `wop/visualization.py`: mesh reconstruction + plotting from planes.
- `wop/cli.py`: command-line runner.
- `tests/`: mandatory correctness tests.
- `docs/walk_on_planes_dirichlet_R3.md`: derivation and usage notes.
