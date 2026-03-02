#pragma once

#include <functional>
#include <string>

namespace wop::estimation {

struct TrajectoryResult {
    double value = 0.0;
    int steps = 0;
    std::string status;  // hit_face | timeout | escaped
};

struct EstimateResult {
    double J = 0.0;
    double S2 = 0.0;
    double eps = 0.0;
    int n_total = 0;
    int n_truncated = 0;
    double mean_steps = 0.0;
};

EstimateResult estimate_from_trajectories(int n_paths, const std::function<TrajectoryResult()>& trace_once);

}  // namespace wop::estimation
