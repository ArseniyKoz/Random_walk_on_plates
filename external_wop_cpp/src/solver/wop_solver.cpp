#include "wop/solver/wop_solver.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "wop/sampling/sampling.hpp"
#include "wop/solver/wop_solver_internal.hpp"

namespace wop::solver {

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
    std::optional<double> r_max,
    RMaxMode r_max_mode,
    double r_max_factor) {
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
    const auto& nu = poly.nu();
    const auto& b = poly.b();
    std::size_t face_idx = 0;
    std::vector<double> d_buffer;
    double eta = 1.0;
    std::optional<double> r_max_sq = std::nullopt;
    if (r_max_mode == RMaxMode::Escape && r_max.has_value()) {
        r_max_sq = (*r_max) * (*r_max);
    }
    const auto r_max_proj = detail::resolve_r_max_projection(poly, x0, r_max, r_max_mode, r_max_factor);
    const auto weighted_boundary_value = [&](double boundary_value) {
        return u_inf + eta * (boundary_value - u_inf);
    };
    try {
        face_idx = poly.closest_outside_face_index(x, eps_in);
    } catch (const std::invalid_argument&) {
        if (poly.is_inside_or_on(x, eps_in)) {
            poly.signed_distances_inplace(x, d_buffer);
            const auto scan = detail::scan_distances(d_buffer, eps_in, 0, eps_plane);
            const std::size_t i0 = scan.argmin_abs;
            return estimation::TrajectoryResult{
                weighted_boundary_value(boundary_f(x, static_cast<int>(i0))),
                0,
                estimation::TrajectoryStatus::HitFace};
        }
        throw std::invalid_argument("x0 must belong to the exterior domain.");
    }

    for (int step = 1; step <= max_steps; ++step) {
        if (!math::is_finite(x) || !std::isfinite(eta)) {
            return estimation::TrajectoryResult{u_inf, step - 1, estimation::TrajectoryStatus::Escaped};
        }

        if (r_max_proj.enabled) {
            const double r = math::norm(x - r_max_proj.center);
            if (!std::isfinite(r)) {
                return estimation::TrajectoryResult{u_inf, step - 1, estimation::TrajectoryStatus::Escaped};
            }
            if (r > r_max_proj.rho1) {
                try {
                    x = detail::sample_far_sphere_step(x, r_max_proj.center, r_max_proj.rho, rng);
                } catch (const std::runtime_error&) {
                    return estimation::TrajectoryResult{u_inf, step, estimation::TrajectoryStatus::Timeout};
                }
                eta *= (r_max_proj.rho / r);
                if (!math::is_finite(x) || !std::isfinite(eta)) {
                    return estimation::TrajectoryResult{u_inf, step, estimation::TrajectoryStatus::Escaped};
                }
                try {
                    face_idx = poly.closest_outside_face_index(x, eps_in);
                } catch (const std::invalid_argument&) {
                    if (poly.is_inside_or_on(x, eps_in)) {
                        poly.signed_distances_inplace(x, d_buffer);
                        const auto scan = detail::scan_distances(d_buffer, eps_in, 0, eps_plane);
                        return estimation::TrajectoryResult{
                            weighted_boundary_value(boundary_f(x, static_cast<int>(scan.argmin_abs))),
                            step,
                            estimation::TrajectoryStatus::HitFace};
                    }
                    throw std::invalid_argument("Far-sphere sample must belong to the exterior domain.");
                }
                continue;
            }
        } else if (r_max_sq.has_value() && math::norm2(x) >= *r_max_sq) {
            return estimation::TrajectoryResult{u_inf, step - 1, estimation::TrajectoryStatus::Escaped};
        }

        math::Vec3 nu_i = nu[face_idx];
        double b_i = b[face_idx];
        double d_i = math::dot(nu_i, x) - b_i;

        if (d_i <= eps_in) {
            poly.signed_distances_inplace(x, d_buffer);
            const auto scan_now = detail::scan_distances(d_buffer, eps_in, face_idx, eps_plane);
            if (!scan_now.any_outside) {
                const std::size_t i_hit = scan_now.argmin_abs;
                return estimation::TrajectoryResult{
                    weighted_boundary_value(boundary_f(x, static_cast<int>(i_hit))),
                    step - 1,
                    estimation::TrajectoryStatus::HitFace};
            }
            face_idx = scan_now.argmin_outside;
            nu_i = nu[face_idx];
            b_i = b[face_idx];
            d_i = math::dot(nu_i, x) - b_i;
        }

        math::Vec3 w;
        try {
            w = sampling::sample_tangent_direction(nu_i, rng, min_abs_denom);
        } catch (const std::runtime_error&) {
            return estimation::TrajectoryResult{u_inf, step, estimation::TrajectoryStatus::Timeout};
        }

        const math::Vec3 p = x - d_i * nu_i;
        const double rho = detail::sample_plane_radius(d_i, rng, min_abs_denom);
        math::Vec3 y = p + rho * w;
        y = y - (math::dot(nu_i, y) - b_i) * nu_i;
        if (!math::is_finite(y)) {
            return estimation::TrajectoryResult{u_inf, step, estimation::TrajectoryStatus::Escaped};
        }

        poly.signed_distances_inplace(y, d_buffer);
        const auto scan = detail::scan_distances(d_buffer, eps_in, face_idx, eps_plane);
        if (scan.hit_target_face) {
            return estimation::TrajectoryResult{
                weighted_boundary_value(boundary_f(y, static_cast<int>(face_idx))),
                step,
                estimation::TrajectoryStatus::HitFace};
        }

        if (!scan.any_outside) {
            const std::size_t i_hit = scan.argmin_abs;
            return estimation::TrajectoryResult{
                weighted_boundary_value(boundary_f(y, static_cast<int>(i_hit))),
                step,
                estimation::TrajectoryStatus::HitFace};
        }

        face_idx = scan.argmin_outside;
        x = y;
    }

    return estimation::TrajectoryResult{u_inf, max_steps, estimation::TrajectoryStatus::Timeout};
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
    std::optional<double> r_max,
    RMaxMode r_max_mode,
    double r_max_factor) {
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
            r_max,
            r_max_mode,
            r_max_factor);
    });
}

}  // namespace wop::solver
