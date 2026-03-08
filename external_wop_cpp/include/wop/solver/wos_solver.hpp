#pragma once

#include <functional>
#include <optional>

#include "wop/estimation/estimation.hpp"
#include "wop/geometry/polyhedron.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"

namespace wop::solver {

using WosBoundaryFunc = std::function<double(const math::Vec3&, std::optional<int>)>;

estimation::TrajectoryResult trace_wos_trajectory(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    double rho_scale = 1.0,
    double rho1_scale = 2.0,
    int max_steps = 1'000'000,
    double u_inf = 0.0);

estimation::EstimateResult estimate_wos(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double delta,
    double rho_scale = 1.0,
    double rho1_scale = 2.0,
    int max_steps = 1'000'000,
    double u_inf = 0.0);

}  // namespace wop::solver

