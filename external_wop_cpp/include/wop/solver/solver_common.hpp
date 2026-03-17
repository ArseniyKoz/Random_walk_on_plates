#pragma once

#include <utility>
#include <vector>

#include "wop/geometry/polyhedron.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"

namespace wop::solver::detail {

double det3_rows(const math::Vec3& r0, const math::Vec3& r1, const math::Vec3& r2);

bool solve_plane_triplet(
    const math::Vec3& n0,
    double b0,
    const math::Vec3& n1,
    double b1,
    const math::Vec3& n2,
    double b2,
    double linear_tol,
    math::Vec3& out_vertex);

std::vector<math::Vec3> compute_polyhedron_vertices(
    const geometry::Polyhedron& poly,
    double linear_tol = 1e-12,
    double inside_tol = 1e-9,
    double dedup_tol = 1e-8);

std::pair<math::Vec3, double> compute_polyhedron_bounding_sphere(
    const geometry::Polyhedron& poly);

math::Vec3 sample_unit_orthogonal(
    const math::Vec3& e,
    rng::Rng& rng,
    double min_norm = 1e-14);

math::Vec3 sample_far_sphere_step(
    const math::Vec3& x,
    const math::Vec3& center,
    double rho,
    rng::Rng& rng);

}  // namespace wop::solver::detail
