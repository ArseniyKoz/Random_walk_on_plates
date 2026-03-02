#pragma once

#include <functional>
#include <optional>

#include "wop/estimation/estimation.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"

namespace wop::solver {

using WosBoundaryFunc = std::function<double(const math::Vec3&, std::optional<int>)>;

estimation::TrajectoryResult trace_wos_box_trajectory(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double eps,
    int max_steps = 1'000'000,
    double u_inf = 0.0,
    std::optional<double> r_max = std::nullopt);

estimation::EstimateResult estimate_wos_box(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double eps,
    int max_steps = 1'000'000,
    double u_inf = 0.0,
    std::optional<double> r_max = std::nullopt);

}  // namespace wop::solver
