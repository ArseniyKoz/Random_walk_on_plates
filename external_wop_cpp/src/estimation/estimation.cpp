#include "wop/estimation/estimation.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace wop::estimation {

namespace {

class EstimateAccumulator {
public:
    void add(const TrajectoryResult& tr) {
        ++count_;
        const double delta = tr.value - mean_;
        mean_ += delta / count_;
        m2_ += delta * (tr.value - mean_);
        steps_sum_ += static_cast<std::int64_t>(tr.steps);

        if (tr.status == TrajectoryStatus::Timeout) {
            ++n_timeout_;
        } else if (tr.status == TrajectoryStatus::Escaped) {
            ++n_escaped_;
        }
    }

    EstimateResult finalize(int n_total) const {
        if (n_total <= 0) {
            throw std::invalid_argument("n_total must be positive.");
        }

        const double j = mean_;
        const double s2 = std::max(m2_ / static_cast<double>(n_total), 0.0);
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
    int count_ = 0;
    double mean_ = 0.0;
    double m2_ = 0.0;
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
