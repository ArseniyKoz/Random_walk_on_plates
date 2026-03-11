#include "wop/config/run_config.hpp"

#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>
#include <vector>

#include "wop/geometry/polyhedron.hpp"
#include "wop/rng/rng.hpp"
#include "wop/solver/wop_solver.hpp"
#include "wop/solver/wos_solver.hpp"

namespace wop::config {

namespace {

double evaluate_function(const FunctionConfig& fn, const math::Vec3& point) {
    switch (fn.kind) {
        case FunctionKind::Constant:
            return fn.value;
        case FunctionKind::X:
            return point.x;
        case FunctionKind::Y:
            return point.y;
        case FunctionKind::Z:
            return point.z;
        case FunctionKind::Coulomb:
            return 1.0 / math::norm(point - fn.source);
    }
    throw std::invalid_argument("Unsupported builtin function kind.");
}

geometry::Polyhedron build_polyhedron(const RuntimeConfig& config) {
    std::vector<geometry::Plane> planes;
    planes.reserve(config.geometry.planes.size());
    for (const auto& plane : config.geometry.planes) {
        planes.push_back(geometry::Plane{plane.p, plane.nu});
    }
    auto poly = geometry::build_polyhedron_from_planes(planes);
    return geometry::orient_normals(poly, config.geometry.interior_point);
}

}  // namespace

const char* method_to_string(Method method) noexcept {
    return (method == Method::Wop) ? "wop" : "wos";
}

ConfigRunResult run_config(const RuntimeConfig& config) {
    ConfigRunResult result;
    result.method = config.method;
    result.x0 = config.x0;

    const geometry::Polyhedron poly = build_polyhedron(config);
    const auto boundary_f = [&](const math::Vec3& y, std::optional<int>) {
        return evaluate_function(config.boundary, y);
    };

    rng::Rng rng(config.seed);
    if (config.method == Method::Wop) {
        result.estimate = solver::estimate_wop(
            poly,
            config.x0,
            boundary_f,
            config.n,
            rng,
            std::nullopt,
            std::nullopt,
            1e-14,
            config.max_steps,
            config.u_inf,
            config.wop->r_max,
            config.wop->r_max_mode,
            config.wop->r_max_factor);
    } else {
        result.estimate = solver::estimate_wos(
            poly,
            config.x0,
            boundary_f,
            config.n,
            rng,
            config.wos->delta,
            config.wos->rho_scale,
            config.wos->rho1_scale,
            config.max_steps,
            config.u_inf);
    }

    if (config.reference.has_value()) {
        const double exact = evaluate_function(*config.reference, config.x0);
        result.exact = exact;
        result.abs_error = std::abs(result.estimate.J - exact);
    }

    return result;
}

std::string format_json_result(const ConfigRunResult& result) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "{"
        << "\"method\":\"" << method_to_string(result.method) << "\","
        << "\"x0\":[" << result.x0.x << "," << result.x0.y << "," << result.x0.z << "],"
        << "\"n_total\":" << result.estimate.n_total << ","
        << "\"n_truncated\":" << result.estimate.n_truncated << ","
        << "\"J\":" << result.estimate.J << ","
        << "\"S2\":" << result.estimate.S2 << ","
        << "\"eps\":" << result.estimate.eps << ","
        << "\"mean_steps\":" << result.estimate.mean_steps;
    if (result.exact.has_value()) {
        out << ",\"exact\":" << *result.exact
            << ",\"abs_error\":" << *result.abs_error;
    }
    out << "}\n";
    return out.str();
}

std::string format_text_result(const ConfigRunResult& result) {
    std::ostringstream out;
    out << "method: " << method_to_string(result.method) << "\n";
    out << "x0: [" << result.x0.x << ", " << result.x0.y << ", " << result.x0.z << "]\n";
    out << "n_total: " << result.estimate.n_total << "\n";
    out << "n_truncated: " << result.estimate.n_truncated << "\n";
    out << std::fixed << std::setprecision(10);
    out << "J: " << result.estimate.J << "\n";
    if (result.exact.has_value()) {
        out << "exact: " << *result.exact << "\n";
        out << "abs_error: " << *result.abs_error << "\n";
    }
    out << "S2: " << result.estimate.S2 << "\n";
    out << "eps: " << result.estimate.eps << "\n";
    out << std::setprecision(2) << "mean_steps: " << result.estimate.mean_steps << "\n";
    return out.str();
}

}  // namespace wop::config
