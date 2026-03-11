#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "wop/config/config.hpp"
#include "wop/config/expression.hpp"
#include "wop/config/run_config.hpp"
#include "wop/geometry/polyhedron.hpp"
#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"
#include "wop/solver/wop_solver.hpp"
#include "wop/solver/wos_solver.hpp"

namespace {

void require(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

bool close(double a, double b, double atol = 1e-12) {
    return std::abs(a - b) <= atol;
}

std::filesystem::path write_temp_config(const std::string& stem, const std::string& body) {
    const auto dir = std::filesystem::temp_directory_path() / "wop_cpp_config_tests";
    std::filesystem::create_directories(dir);
    const auto path = dir / (stem + ".yaml");
    std::ofstream out(path, std::ios::binary);
    out << body;
    out.close();
    return path;
}

std::string read_command_output(const std::string& command) {
#if defined(_WIN32)
    FILE* pipe = _popen(command.c_str(), "r");
#else
    FILE* pipe = popen(command.c_str(), "r");
#endif
    if (pipe == nullptr) {
        throw std::runtime_error("Could not launch subprocess.");
    }

    std::string output;
    char buffer[256];
    while (std::fgets(buffer, static_cast<int>(sizeof(buffer)), pipe) != nullptr) {
        output += buffer;
    }

#if defined(_WIN32)
    const int rc = _pclose(pipe);
#else
    const int rc = pclose(pipe);
#endif
    if (rc != 0) {
        throw std::runtime_error("Subprocess returned non-zero exit code.");
    }
    return output;
}

std::filesystem::path find_cli_binary(const char* argv0) {
    const auto self = std::filesystem::absolute(argv0);
    const auto dir = self.parent_path();
#if defined(_WIN32)
    const auto candidate = dir / "wop_cli.exe";
#else
    const auto candidate = dir / "wop_cli";
#endif
    if (!std::filesystem::exists(candidate)) {
        throw std::runtime_error("Could not locate wop_cli next to test executable.");
    }
    return candidate;
}

std::string cube_config_yaml(const std::string& method, bool with_reference) {
    std::string yaml =
        "command: wop_cli --config examples/" + method + "_cube.yaml --json\n"
        "method: " + method + "\n"
        "x0: [3.0, 0.0, 0.0]\n"
        "n: 300\n"
        "seed: 12345\n"
        "max_steps: 200000\n"
        "u_inf: 0.0\n"
        "geometry:\n"
        "  interior_point: [0.0, 0.0, 0.0]\n"
        "  planes:\n"
        "    - p: [1.0, 0.0, 0.0]\n"
        "      nu: [1.0, 0.0, 0.0]\n"
        "    - p: [-1.0, 0.0, 0.0]\n"
        "      nu: [-1.0, 0.0, 0.0]\n"
        "    - p: [0.0, 1.0, 0.0]\n"
        "      nu: [0.0, 1.0, 0.0]\n"
        "    - p: [0.0, -1.0, 0.0]\n"
        "      nu: [0.0, -1.0, 0.0]\n"
        "    - p: [0.0, 0.0, 1.0]\n"
        "      nu: [0.0, 0.0, 1.0]\n"
        "    - p: [0.0, 0.0, -1.0]\n"
        "      nu: [0.0, 0.0, -1.0]\n"
        "boundary:\n"
        "  expression: \"1 / sqrt((x - 0.2)^2 + (y + 0.1)^2 + (z - 0.3)^2)\"\n";
    if (with_reference) {
        yaml +=
            "reference:\n"
            "  expression: \"1 / sqrt((x - 0.2)^2 + (y + 0.1)^2 + (z - 0.3)^2)\"\n";
    }
    if (method == "wop") {
        yaml +=
            "wop:\n"
            "  r_max: 1000000.0\n"
            "  r_max_mode: escape\n"
            "  r_max_factor: 3.0\n";
    } else {
        yaml +=
            "wos:\n"
            "  delta: 1e-3\n"
            "  rho_scale: 1.0\n"
            "  rho1_scale: 2.0\n";
    }
    return yaml;
}

void test_expression_parser_handles_basic_operators_and_functions() {
    const auto expr = wop::config::compile_expression("sqrt((x - 1)^2) + abs(y) + cos(z) + pi - e");
    const double value = expr.evaluate(wop::math::Vec3{4.0, -2.0, 0.0});
    const double expected = 3.0 + 2.0 + 1.0 + std::acos(-1.0) - std::exp(1.0);
    require(close(value, expected, 1e-12), "expression evaluation mismatch");
}

void test_expression_parser_rejects_unknown_identifier() {
    bool caught = false;
    try {
        static_cast<void>(wop::config::compile_expression("x + face"));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "expression parser should reject unknown identifiers");
}

void test_yaml_config_loader_reads_general_plane_geometry() {
    const auto path = write_temp_config("loader_general_poly", cube_config_yaml("wop", true));
    const auto cfg = wop::config::load_config_file(path);

    require(cfg.method == wop::config::Method::Wop, "method mismatch");
    require(cfg.geometry.planes.size() == 6, "plane count mismatch");
    require(close(cfg.x0.x, 3.0) && close(cfg.x0.y, 0.0) && close(cfg.x0.z, 0.0), "x0 mismatch");
    require(cfg.reference_expression.has_value(), "reference expression should be present");
    require(cfg.wop.has_value(), "wop section should be present");
    require(!cfg.wos.has_value(), "wos section should not be present");
}

void test_yaml_config_loader_allows_optional_command_metadata() {
    const auto path = write_temp_config("loader_command_metadata", cube_config_yaml("wop", false));
    const auto cfg = wop::config::load_config_file(path);

    require(cfg.method == wop::config::Method::Wop, "method mismatch for config with command metadata");
    require(cfg.wop.has_value(), "wop section should be present for config with command metadata");
}

void test_yaml_config_loader_rejects_missing_interior_point() {
    const auto path = write_temp_config(
        "loader_missing_interior",
        "method: wop\n"
        "x0: [3.0, 0.0, 0.0]\n"
        "n: 10\n"
        "seed: 1\n"
        "max_steps: 1000\n"
        "geometry:\n"
        "  planes:\n"
        "    - p: [1.0, 0.0, 0.0]\n"
        "      nu: [1.0, 0.0, 0.0]\n"
        "boundary:\n"
        "  expression: \"x\"\n"
        "wop:\n"
        "  r_max: 10.0\n"
        "  r_max_mode: escape\n"
        "  r_max_factor: 3.0\n");

    bool caught = false;
    try {
        static_cast<void>(wop::config::load_config_file(path));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "loader should reject missing interior_point");
}

void test_run_config_wop_matches_direct_solver_statistics() {
    const auto path = write_temp_config("runner_wop", cube_config_yaml("wop", true));
    const auto cfg = wop::config::load_config_file(path);
    const auto run = wop::config::run_config(cfg);

    std::vector<wop::geometry::Plane> planes;
    planes.reserve(cfg.geometry.planes.size());
    for (const auto& plane : cfg.geometry.planes) {
        planes.push_back(wop::geometry::Plane{plane.p, plane.nu});
    }
    auto poly = wop::geometry::build_polyhedron_from_planes(planes);
    poly = wop::geometry::orient_normals(poly, cfg.geometry.interior_point);

    const auto exact_u = [](const wop::math::Vec3& x) {
        return 1.0 / wop::math::norm(x - wop::math::Vec3{0.2, -0.1, 0.3});
    };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };

    wop::rng::Rng rng(cfg.seed);
    const auto direct = wop::solver::estimate_wop(
        poly,
        cfg.x0,
        boundary_f,
        cfg.n,
        rng,
        std::nullopt,
        std::nullopt,
        1e-14,
        cfg.max_steps,
        cfg.u_inf,
        cfg.wop->r_max,
        cfg.wop->r_max_mode,
        cfg.wop->r_max_factor);

    require(close(run.estimate.J, direct.J, 1e-12), "config WOP J mismatch");
    require(run.exact.has_value(), "config WOP should provide exact");
    require(run.abs_error.has_value(), "config WOP should provide abs_error");
}

void test_run_config_wos_matches_direct_solver_statistics_without_reference() {
    const auto path = write_temp_config("runner_wos", cube_config_yaml("wos", false));
    const auto cfg = wop::config::load_config_file(path);
    const auto run = wop::config::run_config(cfg);

    std::vector<wop::geometry::Plane> planes;
    planes.reserve(cfg.geometry.planes.size());
    for (const auto& plane : cfg.geometry.planes) {
        planes.push_back(wop::geometry::Plane{plane.p, plane.nu});
    }
    auto poly = wop::geometry::build_polyhedron_from_planes(planes);
    poly = wop::geometry::orient_normals(poly, cfg.geometry.interior_point);

    const auto exact_u = [](const wop::math::Vec3& x) {
        return 1.0 / wop::math::norm(x - wop::math::Vec3{0.2, -0.1, 0.3});
    };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };

    wop::rng::Rng rng(cfg.seed);
    const auto direct = wop::solver::estimate_wos(
        poly,
        cfg.x0,
        boundary_f,
        cfg.n,
        rng,
        cfg.wos->delta,
        cfg.wos->rho_scale,
        cfg.wos->rho1_scale,
        cfg.max_steps,
        cfg.u_inf);

    require(close(run.estimate.J, direct.J, 1e-12), "config WoS J mismatch");
    require(!run.exact.has_value(), "config WoS should omit exact when reference is absent");
    require(!run.abs_error.has_value(), "config WoS should omit abs_error when reference is absent");
}

void test_cli_accepts_config_and_emits_json(const char* argv0) {
    const auto cfg_path = write_temp_config("cli_wop", cube_config_yaml("wop", true));
    const auto cli_path = find_cli_binary(argv0);
    const std::string cmd =
        "cmd /C \"\"" + cli_path.string() + "\" --config \"" + cfg_path.string() + "\" --json\"";
    const std::string out = read_command_output(cmd);

    require(out.find("\"method\":\"wop\"") != std::string::npos, "CLI JSON should contain method");
    require(out.find("\"exact\":") != std::string::npos, "CLI JSON should contain exact");
    require(out.find("\"abs_error\":") != std::string::npos, "CLI JSON should contain abs_error");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        test_expression_parser_handles_basic_operators_and_functions();
        test_expression_parser_rejects_unknown_identifier();
        test_yaml_config_loader_reads_general_plane_geometry();
        test_yaml_config_loader_allows_optional_command_metadata();
        test_yaml_config_loader_rejects_missing_interior_point();
        test_run_config_wop_matches_direct_solver_statistics();
        test_run_config_wos_matches_direct_solver_statistics_without_reference();
        test_cli_accepts_config_and_emits_json(argc > 0 ? argv[0] : "");
        std::cout << "wop_config_runtime_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_config_runtime_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
