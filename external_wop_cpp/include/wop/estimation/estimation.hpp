#pragma once

#include <functional>

namespace wop::estimation {

enum class TrajectoryStatus { HitFace, Timeout, Escaped };

inline const char* to_string(TrajectoryStatus s) noexcept {
    switch (s) {
        case TrajectoryStatus::HitFace:  return "hit_face";
        case TrajectoryStatus::Timeout:  return "timeout";
        case TrajectoryStatus::Escaped:  return "escaped";
    }
    return "unknown";
}

struct TrajectoryResult {
    double value = 0.0;
    int steps = 0;
    TrajectoryStatus status = TrajectoryStatus::HitFace;
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
