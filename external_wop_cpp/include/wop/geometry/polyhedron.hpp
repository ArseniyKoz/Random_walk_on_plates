#pragma once

#include <cstddef>
#include <vector>

#include "wop/math/vec3.hpp"

namespace wop::geometry {

struct Plane {
    math::Vec3 p;
    math::Vec3 nu;
};

class Polyhedron {
public:
    Polyhedron(std::vector<math::Vec3> nu, std::vector<double> b);

    const std::vector<math::Vec3>& nu() const noexcept { return nu_; }
    const std::vector<double>& b() const noexcept { return b_; }

    std::size_t num_faces() const noexcept;
    void signed_distances_inplace(const math::Vec3& x, std::vector<double>& out) const;
    std::vector<double> signed_distances(const math::Vec3& x) const;
    bool is_inside_or_on(const math::Vec3& x, double eps_in = 0.0) const;
    std::size_t closest_outside_face_index(const math::Vec3& x, double eps_in = 0.0) const;
    double characteristic_length() const noexcept;

private:
    std::vector<math::Vec3> nu_;
    std::vector<double> b_;
};

Polyhedron build_polyhedron_from_planes(const std::vector<Plane>& planes);
Polyhedron orient_normals(const Polyhedron& poly, const math::Vec3& interior_point);

}  // namespace wop::geometry
