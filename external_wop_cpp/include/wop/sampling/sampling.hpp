#pragma once

#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"

namespace wop::sampling {

struct PlaneHitSample {
    math::Vec3 y;
    math::Vec3 omega;
    double t = 0.0;
};

math::Vec3 sample_unit_sphere(rng::Rng& rng);
math::Vec3 sample_tangent_direction(const math::Vec3& nu, rng::Rng& rng, double min_norm = 1e-14, int max_resample = 10000);
PlaneHitSample sample_hit_on_plane_from_point(
    const math::Vec3& x,
    const math::Vec3& nu,
    double b,
    rng::Rng& rng,
    double min_abs_denom = 1e-14,
    int max_resample = 10000);

}  // namespace wop::sampling
