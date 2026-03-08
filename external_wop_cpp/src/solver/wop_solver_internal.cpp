#include "wop/solver/wop_solver_internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "wop/sampling/sampling.hpp"
#include "wop/solver/wop_solver.hpp"

namespace wop::solver::detail {

namespace {

constexpr double kMinOrthNorm = 1e-14;

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
    double linear_tol,
    math::Vec3& out_vertex) {
    const double det_a = det3_rows(n0, n1, n2);
    if (std::abs(det_a) <= linear_tol) {
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

std::vector<math::Vec3> compute_polyhedron_vertices(
    const geometry::Polyhedron& poly,
    double linear_tol = 1e-12,
    double inside_tol = 1e-9,
    double dedup_tol = 1e-8) {
    if (linear_tol <= 0.0 || inside_tol <= 0.0 || dedup_tol <= 0.0) {
        throw std::invalid_argument("All tolerances must be positive.");
    }

    const auto& nu = poly.nu();
    const auto& b = poly.b();
    const std::size_t m = poly.num_faces();
    if (m < 4) {
        return {};
    }

    const double dedup_tol_sq = dedup_tol * dedup_tol;
    std::vector<math::Vec3> vertices;

    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = i + 1; j < m; ++j) {
            for (std::size_t k = j + 1; k < m; ++k) {
                math::Vec3 v{};
                if (!solve_plane_triplet(nu[i], b[i], nu[j], b[j], nu[k], b[k], linear_tol, v)) {
                    continue;
                }

                bool inside = true;
                for (std::size_t face = 0; face < m; ++face) {
                    const double d_face = math::dot(nu[face], v) - b[face];
                    if (d_face > inside_tol) {
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

math::Vec3 sample_unit_orthogonal(const math::Vec3& e, rng::Rng& rng, double min_norm) {
    if (min_norm <= 0.0 || !std::isfinite(min_norm)) {
        throw std::invalid_argument("min_norm must be finite and positive.");
    }
    for (int attempt = 0; attempt < 512; ++attempt) {
        const math::Vec3 omega = sampling::sample_unit_sphere(rng);
        const math::Vec3 w_raw = omega - math::dot(omega, e) * e;
        const double n = math::norm(w_raw);
        if (n > min_norm) {
            return w_raw / n;
        }
    }
    throw std::runtime_error("Could not sample orthogonal direction.");
}

}  // namespace

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
        rho * rho,
        rho1,
    };
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

    const math::Vec3 w = sample_unit_orthogonal(e, rng, kMinOrthNorm);
    const double sin_part = std::sqrt(std::max(1.0 - z1 * z1, 0.0));
    const math::Vec3 dir = z1 * e + sin_part * w;
    return center + rho * dir;
}

}  // namespace wop::solver::detail
