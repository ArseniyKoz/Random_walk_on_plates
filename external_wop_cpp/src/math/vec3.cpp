#include "wop/math/vec3.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wop::math {

Vec3 operator+(const Vec3& a, const Vec3& b) noexcept {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b) noexcept {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator-(const Vec3& v) noexcept {
    return Vec3{-v.x, -v.y, -v.z};
}

Vec3 operator*(const Vec3& v, double s) noexcept {
    return Vec3{v.x * s, v.y * s, v.z * s};
}

Vec3 operator*(double s, const Vec3& v) noexcept {
    return Vec3{v.x * s, v.y * s, v.z * s};
}

Vec3 operator/(const Vec3& v, double s) noexcept {
    return Vec3{v.x / s, v.y / s, v.z / s};
}

double dot(const Vec3& a, const Vec3& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm2(const Vec3& v) noexcept {
    return dot(v, v);
}

double norm(const Vec3& v) noexcept {
    return std::sqrt(norm2(v));
}

bool is_finite(const Vec3& v) noexcept {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

Vec3 normalized(const Vec3& v, double min_norm) {
    if (min_norm < 0.0) {
        throw std::invalid_argument("min_norm must be non-negative.");
    }
    const double n = norm(v);
    if (!std::isfinite(n) || n <= min_norm) {
        throw std::invalid_argument("Cannot normalize vector with tiny or non-finite norm.");
    }
    return v / n;
}

Vec3 clip(const Vec3& v, const Vec3& mn, const Vec3& mx) noexcept {
    return Vec3{
        std::clamp(v.x, mn.x, mx.x),
        std::clamp(v.y, mn.y, mx.y),
        std::clamp(v.z, mn.z, mx.z),
    };
}

double get_component(const Vec3& v, std::size_t idx) {
    switch (idx) {
        case 0:
            return v.x;
        case 1:
            return v.y;
        case 2:
            return v.z;
        default:
            throw std::out_of_range("Vec3 component index out of range.");
    }
}

void set_component(Vec3& v, std::size_t idx, double value) {
    switch (idx) {
        case 0:
            v.x = value;
            break;
        case 1:
            v.y = value;
            break;
        case 2:
            v.z = value;
            break;
        default:
            throw std::out_of_range("Vec3 component index out of range.");
    }
}

}  // namespace wop::math
