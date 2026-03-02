#pragma once

#include <cstddef>
#include <vector>

#include "wop/rng/rng.hpp"

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

}  // namespace wop::solver::detail
