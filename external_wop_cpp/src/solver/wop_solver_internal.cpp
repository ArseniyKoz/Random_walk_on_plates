#include "wop/solver/wop_solver_internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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

}  // namespace wop::solver::detail
