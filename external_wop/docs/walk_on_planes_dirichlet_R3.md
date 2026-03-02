# Walk on Planes for Exterior Dirichlet Problem in R^3

## 1. Problem Statement

Let `Q` be a convex polyhedron in `R^3`, with boundary `Gamma = dQ`, and exterior domain

`D = R^3 \ overline(Q)`.

We want to estimate `u(x0)` for `x0 in D`, where

- `Delta u = 0` in `D`
- `u = f` on `Gamma`
- `u(x) -> u_inf` as `||x|| -> infinity` (typically `u_inf = 0`)

Polyhedron faces are represented only by planes `(p_i, nu_i)`, where `nu_i` is outward unit normal. Internally we store:

- `nu` with shape `(M, 3)`
- `b` with shape `(M,)`, where `b_i = dot(nu_i, p_i)`

Plane equation: `dot(nu_i, x) = b_i`.
Signed distance (for unit normal): `d_i(x) = dot(nu_i, x) - b_i`.
Inside or on boundary: `d_i(x) <= 0` for all `i`.

## 2. WOP Step Kernel

At current exterior point `x`:

1. Active plane selection:
   `i = argmin_{k: d_k(x) > 0} d_k(x)`.
2. Orthogonal projection of `x` to active plane:
   `p = x - d_i(x) * nu_i`.
3. Sample isotropic tangent direction `w` on that plane:
   - draw `g ~ N(0, I3)`
   - project `g` to tangent plane:
     `w_raw = g - (dot(g, nu_i)) * nu_i`
   - normalize: `w = w_raw / ||w_raw||`
   If `||w_raw||` is too small, resample.
4. Sample radial shift `rho` from inverse CDF of Poisson kernel:
   - `alpha ~ Unif(0,1)`
   - `rho = d_i(x) * sqrt(1/(1-alpha)^2 - 1)`
5. Candidate point on active plane:
   `y = p + rho * w`.
6. Numeric stabilization (projection):
   `y := y - (dot(nu_i, y) - b_i) * nu_i`.
7. Face hit test without vertices/edges:
   - `|d_i(y)| <= eps_plane`
   - `d_j(y) <= eps_in` for all `j`
   If true, trajectory stops and returns `f(y)`.
8. Otherwise, continue from `x := y` with next active plane
   `i_next = argmin_{j: d_j(y) > eps_in} d_j(y)`.

Optional practical stops:

- `max_steps` timeout
- `R_max` escape radius, returning `u_inf`

## 3. Why Face Membership Test Works

For convex polyhedron represented as intersection of half-spaces,

`Q = intersect_j { x : d_j(x) <= 0 }`.

Face `F_i` is exactly points in plane `i` that satisfy all inequalities:

`F_i = { x : d_i(x)=0 and d_j(x)<=0 for all j }`.

So `|d_i(y)| <= eps_plane` and `d_j(y)<=eps_in` for all `j` is a tolerance version of `y in F_i`.

## 4. Module Overview

- `wop/geometry.py`
  - `Polyhedron(nu, b)`
  - `build_polyhedron_from_planes(...)`
  - `orient_normals(poly, interior_point)`
  - `make_axis_aligned_box(...)`
  - box distance/projection helpers
- `wop/sampling.py`
  - `sample_unit_sphere(rng)`
  - `sample_tangent_direction(...)`
  - `sample_hit_on_plane_from_point(...)`
- `wop/wop.py`
  - `trace_wop_trajectory(...)`
  - `estimate_wop(...)`
  - result dataclasses (`TrajectoryResult`, `EstimateResult`)
- `wop/wos_box.py`
  - `trace_wos_box_trajectory(...)`
  - `estimate_wos_box(...)`
- `wop/visualization.py`
  - `compute_polyhedron_vertices(...)`
  - `build_face_polygons(...)`
  - `build_polyhedron_mesh(...)`
  - `plot_polyhedron(...)` (requires matplotlib)
- `wop/cli.py`
  - runnable example

## 5. Default Parameters

- `eps_in = 1e-12 * L`
- `eps_plane = 1e-12 * L`
- `min_abs_denom = 1e-14`
- `max_steps = 1_000_000`
- WOS test tolerance `eps = 1e-3 * L`

`L` is a characteristic scale from geometry (`Polyhedron.characteristic_length`).

## 6. How to Run

Install:

```bash
pip install -e .[test]
```

Tests:

```bash
pytest -q
```

Optional visualization dependency:

```bash
pip install -e .[viz]
```

CLI example:

```bash
python -m wop.cli --example box --x0 "3 0 0" --n 50000
```

## 7. Tests Included

1. Plane sampler quantiles for `z=0`, `x=(0,0,h)` against
   `F_R(r) = 1 - h / sqrt(h^2 + r^2)`.
2. Exact harmonic solution outside box with
   `u(x) = 1 / ||x-a||`, `a` inside box.
3. WOP-vs-WOS consistency on same box and start points.
4. Mesh reconstruction for box: expected 8 vertices and 6 quadrilateral faces.
