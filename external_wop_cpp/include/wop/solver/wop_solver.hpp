#pragma once

#include <functional>
#include <optional>

#include "wop/estimation/estimation.hpp"
#include "wop/geometry/polyhedron.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"

namespace wop::solver {

using BoundaryFunc = std::function<double(const math::Vec3&, std::optional<int>)>;

enum class RMaxMode {
    Escape,
    Project,
};

estimation::TrajectoryResult trace_wop_trajectory(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const BoundaryFunc& boundary_f,
    rng::Rng& rng,
    double eps_in,
    double eps_plane,
    double min_abs_denom = 1e-14,
    int max_steps = 1'000'000,
    double u_inf = 0.0,
    std::optional<double> r_max = std::nullopt,
    RMaxMode r_max_mode = RMaxMode::Project,
    double r_max_factor = 3.0);

estimation::EstimateResult estimate_wop(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const BoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    std::optional<double> eps_in = std::nullopt,
    std::optional<double> eps_plane = std::nullopt,
    double min_abs_denom = 1e-14,
    int max_steps = 1'000'000,
    double u_inf = 0.0,
    std::optional<double> r_max = std::nullopt,
    RMaxMode r_max_mode = RMaxMode::Project,
    double r_max_factor = 3.0);

}  // namespace wop::solver
