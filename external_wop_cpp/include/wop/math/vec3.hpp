#pragma once

#include <cstddef>

namespace wop::math {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b) noexcept;
Vec3 operator-(const Vec3& a, const Vec3& b) noexcept;
Vec3 operator-(const Vec3& v) noexcept;
Vec3 operator*(const Vec3& v, double s) noexcept;
Vec3 operator*(double s, const Vec3& v) noexcept;
Vec3 operator/(const Vec3& v, double s) noexcept;

double dot(const Vec3& a, const Vec3& b) noexcept;
double norm2(const Vec3& v) noexcept;
double norm(const Vec3& v) noexcept;
bool is_finite(const Vec3& v) noexcept;
Vec3 normalized(const Vec3& v, double min_norm = 0.0);
Vec3 clip(const Vec3& v, const Vec3& mn, const Vec3& mx) noexcept;

double get_component(const Vec3& v, std::size_t idx);
void set_component(Vec3& v, std::size_t idx, double value);

}  // namespace wop::math
