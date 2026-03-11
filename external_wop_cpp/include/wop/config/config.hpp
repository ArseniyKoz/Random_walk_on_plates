#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

#include "wop/config/expression.hpp"
#include "wop/math/vec3.hpp"
#include "wop/solver/wop_solver.hpp"

namespace wop::config {

enum class Method {
    Wop,
    Wos,
};

struct PlaneConfig {
    math::Vec3 p{};
    math::Vec3 nu{};
};

struct GeometryConfig {
    math::Vec3 interior_point{};
    std::vector<PlaneConfig> planes;
};

struct WopConfig {
    std::optional<double> r_max = std::nullopt;
    solver::RMaxMode r_max_mode = solver::RMaxMode::Escape;
    double r_max_factor = 3.0;
};

struct WosConfig {
    double delta = 1e-3;
    double rho_scale = 1.0;
    double rho1_scale = 2.0;
};

struct RuntimeConfig {
    Method method = Method::Wop;
    math::Vec3 x0{};
    int n = 0;
    std::uint64_t seed = 0;
    int max_steps = 0;
    double u_inf = 0.0;
    GeometryConfig geometry;
    CompiledExpression boundary_expression;
    std::optional<CompiledExpression> reference_expression = std::nullopt;
    std::optional<WopConfig> wop = std::nullopt;
    std::optional<WosConfig> wos = std::nullopt;
};

RuntimeConfig load_config_file(const std::filesystem::path& path);

}  // namespace wop::config
