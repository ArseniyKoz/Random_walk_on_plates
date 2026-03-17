#include "wop/solver/wop_solver_internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "wop/solver/wop_solver.hpp"

namespace wop::solver::detail {

DistanceScanResult scan_distances(
    const std::vector<double>& values,
    double eps_in,
    std::size_t target_idx,
    double eps_plane) {
    if (values.empty()) {
        throw std::invalid_argument("scan_distances requires non-empty input.");
    }
    if (target_idx >= values.size()) {
        throw std::invalid_argument("target_idx is out of range.");
    }

    DistanceScanResult result;
    double best_abs = std::abs(values[0]);
    double best_outside = std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < values.size(); ++i) {
        const double value = values[i];
        if (value > eps_in) {
            result.any_outside = true;
            result.all_inside = false;
            if (value < best_outside) {
                best_outside = value;
                result.argmin_outside = i;
            }
        } else if (!std::isfinite(value)) {
            result.all_inside = false;
        }

        const double abs_value = std::abs(value);
        if (abs_value < best_abs) {
            best_abs = abs_value;
            result.argmin_abs = i;
        }
    }

    result.hit_target_face = (std::abs(values[target_idx]) <= eps_plane) && result.all_inside;
    return result;
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

RMaxProjectionConfig resolve_r_max_projection(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    std::optional<double> r_max,
    solver::RMaxMode r_max_mode,
    double r_max_factor) {
    if (r_max.has_value() && *r_max <= 0.0) {
        throw std::invalid_argument("r_max must be positive when provided.");
    }

    if (r_max_mode == solver::RMaxMode::Escape) {
        return RMaxProjectionConfig{};
    }
    if (r_max_mode != solver::RMaxMode::Project) {
        throw std::invalid_argument("Unsupported r_max_mode.");
    }
    if (r_max_factor <= 1.0 || !std::isfinite(r_max_factor)) {
        throw std::invalid_argument("r_max_factor must be finite and greater than 1.0.");
    }

    const auto [center, poly_radius] = compute_polyhedron_bounding_sphere(poly);
    double rho = r_max.value_or(r_max_factor * poly_radius);
    if (!r_max.has_value()) {
        const double dist_x0 = math::norm(x0 - center);
        const double eps_floor = 1e-12 * std::max(1.0, poly.characteristic_length());
        const double min_rho = std::max(1.01 * dist_x0, eps_floor);
        rho = std::max(rho, min_rho);
    }
    if (!std::isfinite(rho) || rho <= 0.0) {
        throw std::invalid_argument("Effective rho must be positive and finite.");
    }
    const double rho1 = r_max_factor * rho;
    if (!std::isfinite(rho1) || !(rho1 > rho)) {
        throw std::invalid_argument("Effective rho1 must be finite and greater than rho.");
    }

    return RMaxProjectionConfig{
        true,
        center,
        rho,
        rho1,
    };
}

}  // namespace wop::solver::detail
