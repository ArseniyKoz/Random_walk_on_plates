#include "wop/geometry/box.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace wop::geometry {

namespace {

bool is_close(double a, double b, double rtol = 1e-5, double atol = 1e-8) {
    return std::abs(a - b) <= (atol + rtol * std::abs(b));
}

}  // namespace

Polyhedron make_axis_aligned_box(const math::Vec3& min_corner, const math::Vec3& max_corner) {
    if (!(max_corner.x > min_corner.x && max_corner.y > min_corner.y && max_corner.z > min_corner.z)) {
        throw std::invalid_argument("Each component of max_corner must exceed min_corner.");
    }

    const std::vector<Plane> planes = {
        Plane{math::Vec3{min_corner.x, 0.0, 0.0}, math::Vec3{-1.0, 0.0, 0.0}},
        Plane{math::Vec3{max_corner.x, 0.0, 0.0}, math::Vec3{1.0, 0.0, 0.0}},
        Plane{math::Vec3{0.0, min_corner.y, 0.0}, math::Vec3{0.0, -1.0, 0.0}},
        Plane{math::Vec3{0.0, max_corner.y, 0.0}, math::Vec3{0.0, 1.0, 0.0}},
        Plane{math::Vec3{0.0, 0.0, min_corner.z}, math::Vec3{0.0, 0.0, -1.0}},
        Plane{math::Vec3{0.0, 0.0, max_corner.z}, math::Vec3{0.0, 0.0, 1.0}},
    };

    return build_polyhedron_from_planes(planes);
}

double distance_to_box(const math::Vec3& x, const math::Vec3& min_corner, const math::Vec3& max_corner) {
    const math::Vec3 delta{
        std::max(std::max(min_corner.x - x.x, 0.0), x.x - max_corner.x),
        std::max(std::max(min_corner.y - x.y, 0.0), x.y - max_corner.y),
        std::max(std::max(min_corner.z - x.z, 0.0), x.z - max_corner.z),
    };
    return math::norm(delta);
}

math::Vec3 closest_point_on_box_boundary(
    const math::Vec3& x,
    const math::Vec3& min_corner,
    const math::Vec3& max_corner) {
    const math::Vec3 c = math::clip(x, min_corner, max_corner);
    const bool outside = (x.x < min_corner.x || x.x > max_corner.x || x.y < min_corner.y || x.y > max_corner.y ||
                          x.z < min_corner.z || x.z > max_corner.z);
    if (outside) {
        return c;
    }

    const bool on_boundary =
        is_close(c.x, min_corner.x) || is_close(c.x, max_corner.x) || is_close(c.y, min_corner.y) ||
        is_close(c.y, max_corner.y) || is_close(c.z, min_corner.z) || is_close(c.z, max_corner.z);
    if (on_boundary) {
        return c;
    }

    const std::array<double, 6> distances = {
        c.x - min_corner.x,
        c.y - min_corner.y,
        c.z - min_corner.z,
        max_corner.x - c.x,
        max_corner.y - c.y,
        max_corner.z - c.z,
    };

    std::size_t argmin = 0;
    for (std::size_t i = 1; i < distances.size(); ++i) {
        if (distances[i] < distances[argmin]) {
            argmin = i;
        }
    }

    math::Vec3 y = c;
    switch (argmin) {
        case 0:
            y.x = min_corner.x;
            break;
        case 1:
            y.y = min_corner.y;
            break;
        case 2:
            y.z = min_corner.z;
            break;
        case 3:
            y.x = max_corner.x;
            break;
        case 4:
            y.y = max_corner.y;
            break;
        case 5:
            y.z = max_corner.z;
            break;
        default:
            throw std::logic_error("Unexpected argmin index.");
    }
    return y;
}

}  // namespace wop::geometry
