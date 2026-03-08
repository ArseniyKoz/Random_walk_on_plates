#include "wop/solver/wos_box_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "wop/geometry/box.hpp"
#include "wop/sampling/sampling.hpp"

namespace wop::solver {

namespace {

constexpr double kLinearTol = 1e-12;
constexpr double kInsideTol = 1e-9;
constexpr double kDedupTol = 1e-8;
constexpr double kLambdaTol = 1e-10;
constexpr double kFeasTol = 1e-9;
constexpr double kMinOrthNorm = 1e-14;

struct ProjectionResult {
    math::Vec3 point{};
    double distance = 0.0;
};

double det3_rows(const math::Vec3& r0, const math::Vec3& r1, const math::Vec3& r2) {
    return r0.x * (r1.y * r2.z - r1.z * r2.y) - r0.y * (r1.x * r2.z - r1.z * r2.x) +
           r0.z * (r1.x * r2.y - r1.y * r2.x);
}

bool solve_plane_triplet(
    const math::Vec3& n0,
    double b0,
    const math::Vec3& n1,
    double b1,
    const math::Vec3& n2,
    double b2,
    math::Vec3& out_vertex) {
    const double det_a = det3_rows(n0, n1, n2);
    if (std::abs(det_a) <= kLinearTol) {
        return false;
    }

    const math::Vec3 c0{b0, n0.y, n0.z};
    const math::Vec3 c1{b1, n1.y, n1.z};
    const math::Vec3 c2{b2, n2.y, n2.z};
    const double det_x = det3_rows(c0, c1, c2);

    const math::Vec3 d0{n0.x, b0, n0.z};
    const math::Vec3 d1{n1.x, b1, n1.z};
    const math::Vec3 d2{n2.x, b2, n2.z};
    const double det_y = det3_rows(d0, d1, d2);

    const math::Vec3 e0{n0.x, n0.y, b0};
    const math::Vec3 e1{n1.x, n1.y, b1};
    const math::Vec3 e2{n2.x, n2.y, b2};
    const double det_z = det3_rows(e0, e1, e2);

    out_vertex = math::Vec3{det_x / det_a, det_y / det_a, det_z / det_a};
    return math::is_finite(out_vertex);
}

std::vector<math::Vec3> compute_polyhedron_vertices(const geometry::Polyhedron& poly) {
    const auto& nu = poly.nu();
    const auto& b = poly.b();
    const std::size_t m = poly.num_faces();
    if (m < 4) {
        return {};
    }

    const double dedup_tol_sq = kDedupTol * kDedupTol;
    std::vector<math::Vec3> vertices;

    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            for (std::size_t k = j + 1; k < m; ++k) {
                math::Vec3 v{};
                if (!solve_plane_triplet(nu[i], b[i], nu[j], b[j], nu[k], b[k], v)) {
                    continue;
                }

                bool inside = true;
                for (std::size_t face = 0; face < m; ++face) {
                    const double d_face = math::dot(nu[face], v) - b[face];
                    if (d_face > kInsideTol) {
                        inside = false;
                        break;
                    }
                }
                if (!inside) {
                    continue;
                }

                bool duplicate = false;
                for (const auto& u : vertices) {
                    if (math::norm2(v - u) <= dedup_tol_sq) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    vertices.push_back(v);
                }
            }
        }
    }

    return vertices;
}

std::pair<math::Vec3, double> compute_polyhedron_bounding_sphere(const geometry::Polyhedron& poly) {
    const auto vertices = compute_polyhedron_vertices(poly);
    if (vertices.empty()) {
        const double fallback = std::max(1.0, poly.characteristic_length());
        return {math::Vec3{0.0, 0.0, 0.0}, fallback};
    }

    math::Vec3 v_min = vertices.front();
    math::Vec3 v_max = vertices.front();
    for (const auto& v : vertices) {
        v_min.x = std::min(v_min.x, v.x);
        v_min.y = std::min(v_min.y, v.y);
        v_min.z = std::min(v_min.z, v.z);
        v_max.x = std::max(v_max.x, v.x);
        v_max.y = std::max(v_max.y, v.y);
        v_max.z = std::max(v_max.z, v.z);
    }

    const math::Vec3 center = 0.5 * (v_min + v_max);
    double radius = 0.0;
    for (const auto& v : vertices) {
        radius = std::max(radius, math::norm(v - center));
    }

    const double eps_floor = 1e-12 * std::max(1.0, poly.characteristic_length());
    radius = std::max(radius, eps_floor);
    return {center, radius};
}

bool solve_small_linear_system(
    const std::vector<std::vector<double>>& a,
    const std::vector<double>& b,
    std::vector<double>& x) {
    const std::size_t n = a.size();
    if (n == 0 || n > 3 || b.size() != n) {
        return false;
    }
    for (const auto& row : a) {
        if (row.size() != n) {
            return false;
        }
    }

    double aug[3][4] = {};
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        double best = std::abs(aug[pivot][col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double cur = std::abs(aug[row][col]);
            if (cur > best) {
                best = cur;
                pivot = row;
            }
        }
        if (best <= 1e-14) {
            return false;
        }
        if (pivot != col) {
            for (std::size_t j = col; j <= n; ++j) {
                std::swap(aug[col][j], aug[pivot][j]);
            }
        }

        const double pivot_val = aug[col][col];
        for (std::size_t j = col; j <= n; ++j) {
            aug[col][j] /= pivot_val;
        }
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = aug[row][col];
            for (std::size_t j = col; j <= n; ++j) {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    x.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = aug[i][n];
    }
    return true;
}

ProjectionResult project_to_polyhedron_boundary(
    const geometry::Polyhedron& poly,
    const math::Vec3& x) {
    const auto& nu = poly.nu();
    const auto& b = poly.b();
    const std::size_t m = poly.num_faces();

    bool found = false;
    ProjectionResult best{};
    double best_dist_sq = std::numeric_limits<double>::infinity();

    auto try_active_set = [&](const std::vector<std::size_t>& active) {
        const std::size_t k = active.size();
        std::vector<std::vector<double>> gram(k, std::vector<double>(k, 0.0));
        std::vector<double> rhs(k, 0.0);
        for (std::size_t i = 0; i < k; ++i) {
            const std::size_t fi = active[i];
            rhs[i] = math::dot(nu[fi], x) - b[fi];
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t fj = active[j];
                gram[i][j] = math::dot(nu[fi], nu[fj]);
            }
        }

        std::vector<double> lambda;
        if (!solve_small_linear_system(gram, rhs, lambda)) {
            return;
        }
        for (double value : lambda) {
            if (value < -kLambdaTol) {
                return;
            }
        }

        math::Vec3 y = x;
        for (std::size_t i = 0; i < k; ++i) {
            y = y - lambda[i] * nu[active[i]];
        }
        if (!math::is_finite(y)) {
            return;
        }

        for (std::size_t face = 0; face < m; ++face) {
            const double d_face = math::dot(nu[face], y) - b[face];
            if (d_face > kFeasTol) {
                return;
            }
        }

        const double dist_sq = math::norm2(x - y);
        if (dist_sq < best_dist_sq) {
            best_dist_sq = dist_sq;
            best.point = y;
            best.distance = std::sqrt(std::max(dist_sq, 0.0));
            found = true;
        }
    };

    for (std::size_t i = 0; i < m; ++i) {
        try_active_set({i});
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            try_active_set({i, j});
        }
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            for (std::size_t k = j + 1; k < m; ++k) {
                try_active_set({i, j, k});
            }
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

math::Vec3 sample_unit_orthogonal(const math::Vec3& e, rng::Rng& rng) {
    for (int attempt = 0; attempt < 512; ++attempt) {
        const math::Vec3 omega = sampling::sample_unit_sphere(rng);
        const math::Vec3 w_raw = omega - math::dot(omega, e) * e;
        const double n = math::norm(w_raw);
        if (n > kMinOrthNorm) {
            return w_raw / n;
        }
    }
    throw std::runtime_error("Could not sample orthogonal direction.");
}

math::Vec3 sample_far_sphere_step(
    const math::Vec3& x,
    const math::Vec3& center,
    double rho,
    rng::Rng& rng) {
    const math::Vec3 dx = x - center;
    const double r = math::norm(dx);
    if (!(r > rho)) {
        throw std::invalid_argument("Far-sphere step requires |x-center| > rho.");
    }
    const math::Vec3 e = dx / r;

    const double alpha = rng.uniform01();
    const double denom = r - rho + 2.0 * alpha * rho;
    if (denom <= 0.0) {
        throw std::runtime_error("Invalid denominator in far-sphere sampling.");
    }

    const double q = (r * r - rho * rho) / denom;
    double z1 = (r * r + rho * rho - q * q) / (2.0 * r * rho);
    z1 = std::clamp(z1, -1.0, 1.0);

    const math::Vec3 w = sample_unit_orthogonal(e, rng);
    const double sin_part = std::sqrt(std::max(1.0 - z1 * z1, 0.0));
    const math::Vec3 dir = z1 * e + sin_part * w;
    return center + rho * dir;
}

estimation::TrajectoryResult trace_wos_poly_trajectory(
    const geometry::Polyhedron& poly,
    const math::Vec3& x0,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    if (delta <= 0.0) {
        throw std::invalid_argument("delta must be positive.");
    }
    if (rho_scale <= 0.0 || !std::isfinite(rho_scale)) {
        throw std::invalid_argument("rho_scale must be finite and positive.");
    }
    if (rho1_scale <= 1.0 || !std::isfinite(rho1_scale)) {
        throw std::invalid_argument("rho1_scale must be finite and greater than 1.0.");
    }
    if (max_steps <= 0) {
        throw std::invalid_argument("max_steps must be positive.");
    }

    math::Vec3 x = x0;
    try {
        static_cast<void>(poly.closest_outside_face_index(x, 0.0));
    } catch (const std::invalid_argument&) {
        if (poly.is_inside_or_on(x, 0.0)) {
            const auto d = poly.signed_distances(x);
            const std::size_t i_hit = argmin_abs(d);
            return estimation::TrajectoryResult{boundary_f(x, static_cast<int>(i_hit)), 0, "hit_face"};
        }
        throw std::invalid_argument("x0 must belong to the exterior domain.");
    }

    const auto [center, base_rho] = compute_polyhedron_bounding_sphere(poly);
    const double rho = rho_scale * base_rho;
    const double rho1 = rho1_scale * rho;
    if (!(rho > 0.0) || !(rho1 > rho)) {
        throw std::invalid_argument("Invalid rho/rho1 configuration.");
    }

    double eta = 1.0;
    for (int step = 1; step <= max_steps; ++step) {
        if (!math::is_finite(x) || !std::isfinite(eta)) {
            return estimation::TrajectoryResult{u_inf, step - 1, "timeout"};
        }

        const ProjectionResult proj = project_to_polyhedron_boundary(poly, x);
        if (proj.distance <= delta) {
            return estimation::TrajectoryResult{eta * boundary_f(proj.point, std::nullopt), step, "hit_face"};
        }

        const double r = math::norm(x - center);
        if (r <= rho1) {
            const math::Vec3 omega = sampling::sample_unit_sphere(rng);
            x = x + proj.distance * omega;
        } else {
            x = sample_far_sphere_step(x, center, rho, rng);
            eta *= (rho / r);
        }
    }

    return estimation::TrajectoryResult{u_inf, max_steps, "timeout"};
}

}  // namespace

estimation::TrajectoryResult trace_wos_box_trajectory(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    auto poly = geometry::make_axis_aligned_box(box_min, box_max);
    const math::Vec3 center = 0.5 * (box_min + box_max);
    poly = geometry::orient_normals(poly, center);
    return trace_wos_poly_trajectory(poly, x0, boundary_f, rng, delta, rho_scale, rho1_scale, max_steps, u_inf);
}

estimation::EstimateResult estimate_wos_box(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    auto poly = geometry::make_axis_aligned_box(box_min, box_max);
    const math::Vec3 center = 0.5 * (box_min + box_max);
    poly = geometry::orient_normals(poly, center);

    return estimation::estimate_from_trajectories(n_paths, [&]() {
        return trace_wos_poly_trajectory(poly, x0, boundary_f, rng, delta, rho_scale, rho1_scale, max_steps, u_inf);
    });
}

}  // namespace wop::solver
