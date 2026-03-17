#include "wop/solver/wos_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "wop/sampling/sampling.hpp"
#include "wop/solver/solver_common.hpp"

namespace wop::solver {

namespace {

constexpr double kLinearTol = 1e-12;
constexpr double kLambdaTol = 1e-10;
constexpr double kFeasTol = 1e-9;

struct ProjectionResult {
    math::Vec3 point{};
    double distance = 0.0;
};

struct ActiveSetData {
    std::array<std::size_t, 3> faces{};
    int k = 0;
    std::array<std::array<double, 3>, 3> inv_gram{};
};

struct WosProjectionContext {
    const geometry::Polyhedron* poly = nullptr;
    const std::vector<math::Vec3>* nu = nullptr;
    const std::vector<double>* b = nullptr;
    std::vector<ActiveSetData> active_sets;
    math::Vec3 center{};
    double rho = 0.0;
    double rho1 = 0.0;
};

bool invert_small_symmetric_matrix(
    int n,
    const std::array<std::array<double, 3>, 3>& a,
    std::array<std::array<double, 3>, 3>& inv) {
    if (n == 1) {
        if (std::abs(a[0][0]) <= kLinearTol) {
            return false;
        }
        inv = {};
        inv[0][0] = 1.0 / a[0][0];
        return true;
    }

    if (n == 2) {
        const double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if (std::abs(det) <= kLinearTol) {
            return false;
        }
        inv = {};
        inv[0][0] = a[1][1] / det;
        inv[0][1] = -a[0][1] / det;
        inv[1][0] = -a[1][0] / det;
        inv[1][1] = a[0][0] / det;
        return true;
    }

    if (n == 3) {
        std::array<std::array<double, 6>, 3> aug{};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                aug[i][j] = a[i][j];
            }
            aug[i][3 + i] = 1.0;
        }

        for (int col = 0; col < 3; ++col) {
            int pivot = col;
            double best = std::abs(aug[pivot][col]);
            for (int row = col + 1; row < 3; ++row) {
                const double cur = std::abs(aug[row][col]);
                if (cur > best) {
                    best = cur;
                    pivot = row;
                }
            }
            if (best <= kLinearTol) {
                return false;
            }
            if (pivot != col) {
                std::swap(aug[col], aug[pivot]);
            }

            const double pivot_val = aug[col][col];
            for (int j = col; j < 6; ++j) {
                aug[col][j] /= pivot_val;
            }
            for (int row = 0; row < 3; ++row) {
                if (row == col) {
                    continue;
                }
                const double factor = aug[row][col];
                for (int j = col; j < 6; ++j) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        inv = {};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                inv[i][j] = aug[i][3 + j];
            }
        }
        return true;
    }

    return false;
}

std::vector<ActiveSetData> build_active_sets(const geometry::Polyhedron& poly) {
    const auto& nu = poly.nu();
    const std::size_t m = poly.num_faces();
    std::vector<ActiveSetData> sets;
    sets.reserve(m + (m * (m - 1)) / 2 + (m * (m - 1) * (m - 2)) / 6);

    auto try_add = [&](const std::array<std::size_t, 3>& faces, int k) {
        std::array<std::array<double, 3>, 3> gram{};
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                gram[i][j] = math::dot(nu[faces[static_cast<std::size_t>(i)]], nu[faces[static_cast<std::size_t>(j)]]);
            }
        }

        std::array<std::array<double, 3>, 3> inv{};
        if (!invert_small_symmetric_matrix(k, gram, inv)) {
            return;
        }

        ActiveSetData data{};
        data.faces = faces;
        data.k = k;
        data.inv_gram = inv;
        sets.push_back(data);
    };

    for (std::size_t i = 0; i < m; ++i) {
        try_add({i, 0, 0}, 1);
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            try_add({i, j, 0}, 2);
        }
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            for (std::size_t k = j + 1; k < m; ++k) {
                try_add({i, j, k}, 3);
            }
        }
    }

    if (sets.empty()) {
        throw std::runtime_error("Could not build active sets for WoS projection.");
    }
    return sets;
}

WosProjectionContext resolve_wos_projection_context(
    const geometry::Polyhedron& poly,
    double rho_scale,
    double rho1_scale) {
    if (rho_scale <= 0.0 || !std::isfinite(rho_scale)) {
        throw std::invalid_argument("rho_scale must be finite and positive.");
    }
    if (rho1_scale <= 1.0 || !std::isfinite(rho1_scale)) {
        throw std::invalid_argument("rho1_scale must be finite and greater than 1.0.");
    }

    const auto [center, base_rho] = detail::compute_polyhedron_bounding_sphere(poly);
    const double rho = rho_scale * base_rho;
    const double rho1 = rho1_scale * rho;
    if (!(rho > 0.0) || !(rho1 > rho)) {
        throw std::invalid_argument("Invalid rho/rho1 configuration.");
    }

    WosProjectionContext ctx{};
    ctx.poly = &poly;
    ctx.nu = &poly.nu();
    ctx.b = &poly.b();
    ctx.active_sets = build_active_sets(poly);
    ctx.center = center;
    ctx.rho = rho;
    ctx.rho1 = rho1;
    return ctx;
}

ProjectionResult project_to_polyhedron_boundary(const WosProjectionContext& ctx, const math::Vec3& x) {
    const auto& nu = *ctx.nu;
    const auto& b = *ctx.b;
    const std::size_t m = ctx.poly->num_faces();

    bool found = false;
    ProjectionResult best{};
    double best_dist_sq = std::numeric_limits<double>::infinity();

    for (const auto& set : ctx.active_sets) {
        std::array<double, 3> rhs{};
        for (int i = 0; i < set.k; ++i) {
            rhs[static_cast<std::size_t>(i)] =
                math::dot(nu[set.faces[static_cast<std::size_t>(i)]], x) - b[set.faces[static_cast<std::size_t>(i)]];
        }

        std::array<double, 3> lambda{};
        for (int i = 0; i < set.k; ++i) {
            double acc = 0.0;
            for (int j = 0; j < set.k; ++j) {
                acc += set.inv_gram[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
                       rhs[static_cast<std::size_t>(j)];
            }
            lambda[static_cast<std::size_t>(i)] = acc;
        }

        bool lambda_ok = true;
        for (int i = 0; i < set.k; ++i) {
            if (lambda[static_cast<std::size_t>(i)] < -kLambdaTol) {
                lambda_ok = false;
                break;
            }
        }
        if (!lambda_ok) {
            continue;
        }

        math::Vec3 y = x;
        for (int i = 0; i < set.k; ++i) {
            y = y - lambda[static_cast<std::size_t>(i)] * nu[set.faces[static_cast<std::size_t>(i)]];
        }
        if (!math::is_finite(y)) {
            continue;
        }

        bool feasible = true;
        for (std::size_t face = 0; face < m; ++face) {
            const double d_face = math::dot(nu[face], y) - b[face];
            if (d_face > kFeasTol) {
                feasible = false;
                break;
            }
        }
        if (!feasible) {
            continue;
        }

        const double dist_sq = math::norm2(x - y);
        if (dist_sq < best_dist_sq) {
            best_dist_sq = dist_sq;
            best.point = y;
            best.distance = std::sqrt(std::max(dist_sq, 0.0));
            found = true;
        }
    }

    if (!found) {
        throw std::runtime_error("Failed to project point to polyhedron boundary.");
    }
    return best;
}

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

estimation::TrajectoryResult trace_wos_with_context(
    const WosProjectionContext& ctx,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    int max_steps,
    double u_inf) {
    if (delta <= 0.0) {
        throw std::invalid_argument("delta must be positive.");
    }
    if (max_steps <= 0) {
        throw std::invalid_argument("max_steps must be positive.");
    }

    const auto& poly = *ctx.poly;
    math::Vec3 x = x0;

    try {
        static_cast<void>(poly.closest_outside_face_index(x, 0.0));
    } catch (const std::invalid_argument&) {
        if (poly.is_inside_or_on(x, 0.0)) {
            const auto d = poly.signed_distances(x);
            const std::size_t i_hit = argmin_abs(d);
            return estimation::TrajectoryResult{boundary_f(x, static_cast<int>(i_hit)), 0, estimation::TrajectoryStatus::HitFace};
        }
        throw std::invalid_argument("x0 must belong to the exterior domain.");
    }

    double eta = 1.0;
    for (int step = 1; step <= max_steps; ++step) {
        if (!math::is_finite(x) || !std::isfinite(eta)) {
            return estimation::TrajectoryResult{u_inf, step - 1, estimation::TrajectoryStatus::Timeout};
        }

        const ProjectionResult proj = project_to_polyhedron_boundary(ctx, x);
        if (proj.distance <= delta) {
            const double bv = boundary_f(proj.point, std::nullopt);
            return estimation::TrajectoryResult{u_inf + eta * (bv - u_inf), step, estimation::TrajectoryStatus::HitFace};
        }

        const double r = math::norm(x - ctx.center);
        if (r <= ctx.rho1) {
            const math::Vec3 omega = sampling::sample_unit_sphere(rng);
            x = x + proj.distance * omega;
        } else {
            x = detail::sample_far_sphere_step(x, ctx.center, ctx.rho, rng);
            eta *= (ctx.rho / r);
        }
    }

    return estimation::TrajectoryResult{u_inf, max_steps, estimation::TrajectoryStatus::Timeout};
}

}  // namespace

estimation::TrajectoryResult trace_wos_trajectory(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    const auto ctx = resolve_wos_projection_context(poly, rho_scale, rho1_scale);
    return trace_wos_with_context(ctx, x0, boundary_f, rng, delta, max_steps, u_inf);
}

estimation::EstimateResult estimate_wos(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    const auto ctx = resolve_wos_projection_context(poly, rho_scale, rho1_scale);
    return estimation::estimate_from_trajectories(n_paths, [&]() {
        return trace_wos_with_context(ctx, x0, boundary_f, rng, delta, max_steps, u_inf);
    });
}

}  // namespace wop::solver
