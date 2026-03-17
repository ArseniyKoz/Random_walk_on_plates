#include "wop/config/problem_functions.hpp"

#include "wop/math/vec3.hpp"

namespace wop::config {

namespace {

double coulomb_a(const math::Vec3& p) {
    const math::Vec3 a{0.2, -0.1, 0.3};
    return 1.0 ;
}

}  // namespace

double boundary_value(const math::Vec3& y) {
    return coulomb_a(y);
}

bool has_reference_value() noexcept {
    return true;
}

double reference_value(const math::Vec3& x) {
    return coulomb_a(x);
}

}  // namespace wop::config
