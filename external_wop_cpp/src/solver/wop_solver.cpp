#include "wop/solver/wop_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "wop/sampling/sampling.hpp"

namespace wop::solver {

namespace {

std::size_t argmin_abs(const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument("argmin_abs requires non-empty input.");
    }
    std::size_t idx = 0;
    double best = std::abs(values[0]);
    for (std::size_t i = 1; i < values.size(); ++i) {
        const double cur = std::abs(values[i]);
        if (cur < best) {
            best = cur;
            idx = i;
        }
    }
    return idx;
}

std::size_t argmin_outside(const std::vector<double>& values, double eps_in) {
    std::size_t idx = values.size();
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (values[i] > eps_in && values[i] < best) {
            best = values[i];
            idx = i;
        }
    }
    if (idx == values.size()) {
        throw std::invalid_argument("No outside face found.");
    }
    return idx;
}

bool any_outside(const std::vector<double>& values, double eps_in) {
    return std::any_of(values.begin(), values.end(), [eps_in](double d) { return d > eps_in; });
}

double sample_plane_radius(double d, rng::Rng& rng, double min_one_minus_alpha) {
    if (d <= 0.0) {
        throw std::invalid_argument("d must be positive.");
    }
    if (min_one_minus_alpha <= 0.0) {
        throw std::invalid_argument("min_one_minus_alpha must be positive.");
    }

    const double alpha = rng.uniform01();
    const double one_minus_alpha = std::max(1.0 - alpha, min_one_minus_alpha);
    const double rho_sq = d * d * (1.0 / (one_minus_alpha * one_minus_alpha) - 1.0);
    return std::sqrt(std::max(rho_sq, 0.0));
}

}  // namespace

estimation::TrajectoryResult trace_wop_trajectory(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const BoundaryFunc& boundary_f,
    rng::Rng& rng,
    double eps_in,
    double eps_plane,
    double min_abs_denom,
    int max_steps,
    double u_inf,
    std::optional<double> r_max) {
    if (max_steps <= 0) {
        throw std::invalid_argument("max_steps must be positive.");
    }
    if (eps_in < 0.0 || eps_plane < 0.0) {
        throw std::invalid_argument("eps_in and eps_plane must be non-negative.");
    }
    if (min_abs_denom <= 0.0) {
        throw std::invalid_argument("min_abs_denom must be positive.");
    }

    math::Vec3 x = x0;
    std::size_t face_idx = 0;
    try {
        face_idx = poly.closest_outside_face_index(x, eps_in);
    } catch (const std::invalid_argument&) {
        if (poly.is_inside_or_on(x, eps_in)) {
            const auto d = poly.signed_distances(x);
            const std::size_t i0 = argmin_abs(d);
            return estimation::TrajectoryResult{boundary_f(x, static_cast<int>(i0)), 0, "hit_face"};
        }
        throw std::invalid_argument("x0 must belong to the exterior domain.");
    }

    for (int step = 1; step <= max_steps; ++step) {
        if (!math::is_finite(x)) {
            return estimation::TrajectoryResult{u_inf, step - 1, "escaped"};
        }
        if (r_max.has_value() && math::norm(x) >= *r_max) {
            return estimation::TrajectoryResult{u_inf, step - 1, "escaped"};
        }

        math::Vec3 nu_i = poly.nu()[face_idx];
        double b_i = poly.b()[face_idx];
        double d_i = math::dot(nu_i, x) - b_i;

        if (d_i <= eps_in) {
            const auto d_now = poly.signed_distances(x);
            if (!any_outside(d_now, eps_in)) {
                const std::size_t i_hit = argmin_abs(d_now);
                return estimation::TrajectoryResult{boundary_f(x, static_cast<int>(i_hit)), step - 1, "hit_face"};
            }
            face_idx = argmin_outside(d_now, eps_in);
            nu_i = poly.nu()[face_idx];
            b_i = poly.b()[face_idx];
            d_i = math::dot(nu_i, x) - b_i;
        }

        math::Vec3 w;
        try {
            w = sampling::sample_tangent_direction(nu_i, rng, min_abs_denom);
        } catch (const std::runtime_error&) {
            return estimation::TrajectoryResult{u_inf, step, "timeout"};
        }

        const math::Vec3 p = x - d_i * nu_i;
        const double rho = sample_plane_radius(d_i, rng, min_abs_denom);
        math::Vec3 y = p + rho * w;
        y = y - (math::dot(nu_i, y) - b_i) * nu_i;
        if (!math::is_finite(y)) {
            return estimation::TrajectoryResult{u_inf, step, "escaped"};
        }

        const auto d = poly.signed_distances(y);
        const bool hit = (std::abs(d[face_idx]) <= eps_plane) &&
                         std::all_of(d.begin(), d.end(), [eps_in](double value) { return value <= eps_in; });
        if (hit) {
            return estimation::TrajectoryResult{boundary_f(y, static_cast<int>(face_idx)), step, "hit_face"};
        }

        if (!any_outside(d, eps_in)) {
            const std::size_t i_hit = argmin_abs(d);
            return estimation::TrajectoryResult{boundary_f(y, static_cast<int>(i_hit)), step, "hit_face"};
        }

        face_idx = argmin_outside(d, eps_in);
        x = y;
    }

    return estimation::TrajectoryResult{u_inf, max_steps, "timeout"};
}

estimation::EstimateResult estimate_wop(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const BoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    std::optional<double> eps_in,
    std::optional<double> eps_plane,
    double min_abs_denom,
    int max_steps,
    double u_inf,
    std::optional<double> r_max) {
    const double L = poly.characteristic_length();
    const double eps_in_eff = eps_in.value_or(1e-12 * L);
    const double eps_plane_eff = eps_plane.value_or(1e-12 * L);

    return estimation::estimate_from_trajectories(n_paths, [&]() {
        return trace_wop_trajectory(
            poly,
            x0,
            boundary_f,
            rng,
            eps_in_eff,
            eps_plane_eff,
            min_abs_denom,
            max_steps,
            u_inf,
            r_max);
    });
}

}  // namespace wop::solver
