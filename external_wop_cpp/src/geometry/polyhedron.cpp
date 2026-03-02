#include "wop/geometry/polyhedron.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace wop::geometry {

namespace {

constexpr double kNormRtol = 1e-10;
constexpr double kNormAtol = 1e-12;

bool is_unit_norm(double value) {
    return std::abs(value - 1.0) <= (kNormAtol + kNormRtol);
}

}  // namespace

Polyhedron::Polyhedron(std::vector<math::Vec3> nu, std::vector<double> b)
    : nu_(std::move(nu)), b_(std::move(b)) {
    if (nu_.empty()) {
        throw std::invalid_argument("Polyhedron must contain at least one face.");
    }
    if (b_.size() != nu_.size()) {
        throw std::invalid_argument("b must have shape (M,) and match nu.");
    }

    for (const auto& n : nu_) {
        const double n_norm = math::norm(n);
        if (n_norm <= 0.0) {
            throw std::invalid_argument("Each face normal must have non-zero norm.");
        }
        if (!is_unit_norm(n_norm)) {
            throw std::invalid_argument("All face normals must be unit vectors.");
        }
    }
}

std::size_t Polyhedron::num_faces() const noexcept {
    return b_.size();
}

void Polyhedron::signed_distances_inplace(const math::Vec3& x, std::vector<double>& out) const {
    out.resize(num_faces());
    for (std::size_t i = 0; i < num_faces(); ++i) {
        out[i] = math::dot(nu_[i], x) - b_[i];
    }
}

std::vector<double> Polyhedron::signed_distances(const math::Vec3& x) const {
    std::vector<double> d;
    signed_distances_inplace(x, d);
    return d;
}

bool Polyhedron::is_inside_or_on(const math::Vec3& x, double eps_in) const {
    const auto d = signed_distances(x);
    return std::all_of(d.begin(), d.end(), [eps_in](double value) { return value <= eps_in; });
}

std::size_t Polyhedron::closest_outside_face_index(const math::Vec3& x, double eps_in) const {
    const auto d = signed_distances(x);
    std::size_t best_i = num_faces();
    double best_d = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < d.size(); ++i) {
        if (d[i] > eps_in && d[i] < best_d) {
            best_d = d[i];
            best_i = i;
        }
    }
    if (best_i == num_faces()) {
        throw std::invalid_argument("Point is not outside with respect to eps_in.");
    }
    return best_i;
}

double Polyhedron::characteristic_length() const noexcept {
    double scale = 1.0;
    for (double value : b_) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

Polyhedron build_polyhedron_from_planes(const std::vector<Plane>& planes) {
    if (planes.empty()) {
        throw std::invalid_argument("At least one plane is required.");
    }

    std::vector<math::Vec3> nu;
    std::vector<double> b;
    nu.reserve(planes.size());
    b.reserve(planes.size());

    for (const auto& plane : planes) {
        const double n_norm = math::norm(plane.nu);
        if (n_norm == 0.0) {
            throw std::invalid_argument("nu has zero norm.");
        }
        const math::Vec3 nu_unit = plane.nu / n_norm;
        nu.push_back(nu_unit);
        b.push_back(math::dot(nu_unit, plane.p));
    }

    return Polyhedron(std::move(nu), std::move(b));
}

Polyhedron orient_normals(const Polyhedron& poly, const math::Vec3& interior_point) {
    std::vector<math::Vec3> nu = poly.nu();
    std::vector<double> b = poly.b();

    for (std::size_t i = 0; i < nu.size(); ++i) {
        const double d = math::dot(nu[i], interior_point) - b[i];
        if (d > 0.0) {
            nu[i] = -nu[i];
            b[i] = -b[i];
        }
    }
    return Polyhedron(std::move(nu), std::move(b));
}

}  // namespace wop::geometry
