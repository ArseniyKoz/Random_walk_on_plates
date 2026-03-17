#pragma once

#include "wop/math/vec3.hpp"

namespace wop::config {

// Edit the implementations in app/problem_functions.cpp to define
// the active problem used by config-driven CLI runs.
double boundary_value(const math::Vec3& y);
bool has_reference_value() noexcept;
double reference_value(const math::Vec3& x);

}  // namespace wop::config
