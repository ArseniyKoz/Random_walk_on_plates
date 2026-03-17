#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "wop/geometry/polyhedron.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"
#include "wop/solver/solver_common.hpp"

namespace wop::solver {

enum class RMaxMode;

}  // namespace wop::solver

namespace wop::solver::detail {

struct DistanceScanResult {
    bool any_outside = false;
    bool all_inside = true;
    std::size_t argmin_abs = 0;
    std::size_t argmin_outside = 0;
    bool hit_target_face = false;
};

DistanceScanResult scan_distances(
    const std::vector<double>& values,
    double eps_in,
    std::size_t target_idx,
    double eps_plane);

double sample_plane_radius(double d, rng::Rng& rng, double min_one_minus_alpha);

struct RMaxProjectionConfig {
    bool enabled = false;
    math::Vec3 center{};
    double rho = 0.0;
    double rho1 = 0.0;
};

RMaxProjectionConfig resolve_r_max_projection(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    std::optional<double> r_max,
    solver::RMaxMode r_max_mode,
    double r_max_factor);

}  // namespace wop::solver::detail
