#include "wop/estimation/estimation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace wop::estimation {

namespace {

class EstimateAccumulator {
public:
    void add(const TrajectoryResult& tr) {
        const double ksi = tr.value;
        j_sum_ += ksi;
        s2_raw_sum_ += ksi * ksi;
        steps_sum_ += static_cast<std::int64_t>(tr.steps);

        if (tr.status == "timeout") {
            ++n_timeout_;
        } else if (tr.status == "escaped") {
            ++n_escaped_;
        }
    }

    EstimateResult finalize(int n_total) const {
        if (n_total <= 0) {
            throw std::invalid_argument("n_total must be positive.");
        }

        const double j = j_sum_ / static_cast<double>(n_total);
        double s2 = s2_raw_sum_ / static_cast<double>(n_total) - j * j;
        s2 = std::max(s2, 0.0);
        const double eps = 3.0 * std::sqrt(s2 / static_cast<double>(n_total));
        const int n_truncated = n_timeout_ + n_escaped_;

        return EstimateResult{
            j,
            s2,
            eps,
            n_total,
            n_truncated,
            static_cast<double>(steps_sum_) / static_cast<double>(n_total),
        };
    }

private:
    double j_sum_ = 0.0;
    double s2_raw_sum_ = 0.0;
    std::int64_t steps_sum_ = 0;
    int n_timeout_ = 0;
    int n_escaped_ = 0;
};

}  // namespace

EstimateResult estimate_from_trajectories(int n_paths, const std::function<TrajectoryResult()>& trace_once) {
    if (n_paths <= 0) {
        throw std::invalid_argument("n_paths must be positive.");
    }

    EstimateAccumulator acc;
    for (int i = 0; i < n_paths; ++i) {
        acc.add(trace_once());
    }
    return acc.finalize(n_paths);
}

}  // namespace wop::estimation
