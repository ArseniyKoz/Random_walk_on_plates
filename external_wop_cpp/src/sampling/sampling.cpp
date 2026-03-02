#include "wop/sampling/sampling.hpp"

#include <cmath>
#include <stdexcept>

namespace wop::sampling {

math::Vec3 sample_unit_sphere(rng::Rng& rng) {
    while (true) {
        const math::Vec3 w = rng.normal3();
        const double n = math::norm(w);
        if (n > 0.0) {
            return w / n;
        }
    }
}

math::Vec3 sample_tangent_direction(const math::Vec3& nu, rng::Rng& rng, double min_norm, int max_resample) {
    if (min_norm <= 0.0) {
        throw std::invalid_argument("min_norm must be positive.");
    }
    if (max_resample <= 0) {
        throw std::invalid_argument("max_resample must be positive.");
    }

    const double nu_norm = math::norm(nu);
    if (nu_norm == 0.0) {
        throw std::invalid_argument("nu must have non-zero norm.");
    }
    const math::Vec3 nu_unit = nu / nu_norm;

    for (int attempt = 0; attempt < max_resample; ++attempt) {
        const math::Vec3 g = rng.normal3();
        const math::Vec3 w_raw = g - math::dot(g, nu_unit) * nu_unit;
        const double n = math::norm(w_raw);
        if (!std::isfinite(n) || n <= min_norm) {
            continue;
        }
        return w_raw / n;
    }

    throw std::runtime_error("Failed to sample stable tangent direction.");
}

PlaneHitSample sample_hit_on_plane_from_point(
    const math::Vec3& x,
    const math::Vec3& nu,
    double b,
    rng::Rng& rng,
    double min_abs_denom,
    int max_resample) {
    if (min_abs_denom <= 0.0) {
        throw std::invalid_argument("min_abs_denom must be positive.");
    }
    if (max_resample <= 0) {
        throw std::invalid_argument("max_resample must be positive.");
    }

    const double d = math::dot(nu, x) - b;
    if (d <= 0.0) {
        throw std::invalid_argument("x must satisfy dot(nu,x)-b > 0 for plane sampling.");
    }

    for (int attempt = 0; attempt < max_resample; ++attempt) {
        math::Vec3 omega = sample_unit_sphere(rng);
        if (math::dot(nu, omega) > 0.0) {
            omega = -omega;
        }

        const double denom = math::dot(nu, omega);
        if (std::abs(denom) < min_abs_denom) {
            continue;
        }

        const double t = -d / denom;
        if (t <= 0.0 || !std::isfinite(t)) {
            continue;
        }

        const math::Vec3 y = x + t * omega;
        if (!math::is_finite(y)) {
            continue;
        }

        return PlaneHitSample{y, omega, t};
    }

    throw std::runtime_error("Failed to sample stable ray-plane intersection.");
}

}  // namespace wop::sampling
