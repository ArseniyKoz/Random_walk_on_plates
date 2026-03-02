#pragma once

#include "wop/geometry/polyhedron.hpp"

namespace wop::geometry {

Polyhedron make_axis_aligned_box(const math::Vec3& min_corner, const math::Vec3& max_corner);
double distance_to_box(const math::Vec3& x, const math::Vec3& min_corner, const math::Vec3& max_corner);
math::Vec3 closest_point_on_box_boundary(
    const math::Vec3& x,
    const math::Vec3& min_corner,
    const math::Vec3& max_corner);

}  // namespace wop::geometry
