#include "wop/solver/wos_box_solver.hpp"

#include <stdexcept>

#include "wop/geometry/box.hpp"
#include "wop/sampling/sampling.hpp"

namespace wop::solver {

estimation::TrajectoryResult trace_wos_box_trajectory(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double eps,
    int max_steps,
    double u_inf,
    std::optional<double> r_max) {
    if (eps <= 0.0) {
        throw std::invalid_argument("eps must be positive.");
    }
    if (max_steps <= 0) {
        throw std::invalid_argument("max_steps must be positive.");
    }

    math::Vec3 x = x0;
    std::optional<double> r_max_sq = std::nullopt;
    if (r_max.has_value()) {
        r_max_sq = (*r_max) * (*r_max);
    }
    for (int step = 1; step <= max_steps; ++step) {
        if (r_max_sq.has_value() && math::norm2(x) >= *r_max_sq) {
            return estimation::TrajectoryResult{u_inf, step - 1, "escaped"};
        }

        const double dist = geometry::distance_to_box(x, box_min, box_max);
        if (dist <= eps) {
            const math::Vec3 y = geometry::closest_point_on_box_boundary(x, box_min, box_max);
            return estimation::TrajectoryResult{boundary_f(y, std::nullopt), step, "hit_face"};
        }

        const math::Vec3 omega = sampling::sample_unit_sphere(rng);
        x = x + dist * omega;
    }

    return estimation::TrajectoryResult{u_inf, max_steps, "timeout"};
}

estimation::EstimateResult estimate_wos_box(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double eps,
    int max_steps,
    double u_inf,
    std::optional<double> r_max) {
    return estimation::estimate_from_trajectories(n_paths, [&]() {
        return trace_wos_box_trajectory(x0, box_min, box_max, boundary_f, rng, eps, max_steps, u_inf, r_max);
    });
}

}  // namespace wop::solver
