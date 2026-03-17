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
#include "wop/config/problem_functions.hpp"
#include "wop/config/run_config.hpp"
#include "wop/geometry/box.hpp"
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

std::filesystem::path find_examples_dir(const char* argv0) {
    const auto self = std::filesystem::absolute(argv0);
    const auto root = self.parent_path().parent_path();
    const auto examples = root / "examples";
    if (!std::filesystem::exists(examples)) {
        throw std::runtime_error("Could not locate examples directory.");
    }
    return examples;
}

std::string cube_config_yaml(const std::string& method) {
    std::string yaml =
        "command: .\\build\\Release\\wop_cli.exe --config examples/" + method + "_cube.yaml --json\n"
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
        "      nu: [0.0, 0.0, -1.0]\n";
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

void test_problem_functions_are_finite_on_reference_point() {
    const wop::math::Vec3 p{3.0, 0.0, 0.0};
    require(std::isfinite(wop::config::boundary_value(p)), "boundary_value must stay finite at x0");
    if (wop::config::has_reference_value()) {
        require(std::isfinite(wop::config::reference_value(p)), "reference_value must stay finite at x0");
    }
}

std::vector<wop::geometry::Plane> unit_cube_planes() {
    return {
        {wop::math::Vec3{1.0, 0.0, 0.0}, wop::math::Vec3{1.0, 0.0, 0.0}},
        {wop::math::Vec3{-1.0, 0.0, 0.0}, wop::math::Vec3{-1.0, 0.0, 0.0}},
        {wop::math::Vec3{0.0, 1.0, 0.0}, wop::math::Vec3{0.0, 1.0, 0.0}},
        {wop::math::Vec3{0.0, -1.0, 0.0}, wop::math::Vec3{0.0, -1.0, 0.0}},
        {wop::math::Vec3{0.0, 0.0, 1.0}, wop::math::Vec3{0.0, 0.0, 1.0}},
        {wop::math::Vec3{0.0, 0.0, -1.0}, wop::math::Vec3{0.0, 0.0, -1.0}},
    };
}

wop::geometry::Polyhedron make_general_cube_polyhedron() {
    auto poly = wop::geometry::build_polyhedron_from_planes(unit_cube_planes());
    return wop::geometry::orient_normals(poly, wop::math::Vec3{0.0, 0.0, 0.0});
}

wop::geometry::Polyhedron make_axis_box_polyhedron() {
    auto poly = wop::geometry::make_axis_aligned_box(
        wop::math::Vec3{-1.0, -1.0, -1.0},
        wop::math::Vec3{1.0, 1.0, 1.0});
    return wop::geometry::orient_normals(poly, wop::math::Vec3{0.0, 0.0, 0.0});
}

void test_yaml_config_loader_accepts_config_without_boundary_section() {
    const auto path = write_temp_config("loader_no_boundary", cube_config_yaml("wop"));
    const auto cfg = wop::config::load_config_file(path);

    require(cfg.method == wop::config::Method::Wop, "method mismatch");
    require(cfg.geometry.planes.size() == 6, "plane count mismatch");
    require(close(cfg.x0.x, 3.0) && close(cfg.x0.y, 0.0) && close(cfg.x0.z, 0.0), "x0 mismatch");
    require(cfg.wop.has_value(), "wop section should be present");
    require(!cfg.wos.has_value(), "wos section should not be present");
}

void test_yaml_config_loader_rejects_legacy_boundary_section() {
    const auto path = write_temp_config(
        "loader_legacy_boundary",
        "method: wop\n"
        "x0: [3.0, 0.0, 0.0]\n"
        "n: 10\n"
        "seed: 1\n"
        "max_steps: 1000\n"
        "geometry:\n"
        "  interior_point: [0.0, 0.0, 0.0]\n"
        "  planes:\n"
        "    - p: [1.0, 0.0, 0.0]\n"
        "      nu: [1.0, 0.0, 0.0]\n"
        "boundary:\n"
        "  kind: coulomb\n"
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
    require(caught, "loader should reject legacy boundary section");
}

void test_yaml_config_loader_allows_optional_command_metadata() {
    const auto path = write_temp_config("loader_command_metadata", cube_config_yaml("wop"));
    const auto cfg = wop::config::load_config_file(path);

    require(cfg.method == wop::config::Method::Wop, "method mismatch for config with command metadata");
    require(cfg.wop.has_value(), "wop section should be present for config with command metadata");
}

void test_example_configs_are_tuned_for_fast_runs(const char* argv0) {
    const auto examples = find_examples_dir(argv0);
    const auto wop_cfg = wop::config::load_config_file(examples / "box_wop.yaml");
    const auto wos_cfg = wop::config::load_config_file(examples / "box_wos.yaml");
    const auto eq_cfg = wop::config::load_config_file(examples / "box_wop_legacy_equivalent.yaml");

    require(wop_cfg.method == wop::config::Method::Wop, "box_wop example must use WOP");
    require(wos_cfg.method == wop::config::Method::Wos, "box_wos example must use WoS");
    require(eq_cfg.method == wop::config::Method::Wop, "legacy-equivalent example must use WOP");
    require(wop_cfg.n <= 50000, "box_wop example should stay interactive");
    require(wos_cfg.n <= 50000, "box_wos example should stay interactive");
    require(eq_cfg.n == 50000, "legacy-equivalent example should mirror legacy sample count");
    require(eq_cfg.max_steps == 1000000, "legacy-equivalent example should mirror legacy max_steps");
    require(eq_cfg.wop.has_value(), "legacy-equivalent example must define WOP parameters");
    require(eq_cfg.wop->r_max.has_value() && close(*eq_cfg.wop->r_max, 1e6), "legacy-equivalent r_max mismatch");
    require(eq_cfg.wop->r_max_mode == wop::solver::RMaxMode::Escape, "legacy-equivalent mode must match CLI default");
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
    const auto path = write_temp_config("runner_wop", cube_config_yaml("wop"));
    const auto cfg = wop::config::load_config_file(path);
    const auto run = wop::config::run_config(cfg);

    std::vector<wop::geometry::Plane> planes;
    planes.reserve(cfg.geometry.planes.size());
    for (const auto& plane : cfg.geometry.planes) {
        planes.push_back(wop::geometry::Plane{plane.p, plane.nu});
    }
    auto poly = wop::geometry::build_polyhedron_from_planes(planes);
    poly = wop::geometry::orient_normals(poly, cfg.geometry.interior_point);

    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return wop::config::boundary_value(y);
    };

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
    require(wop::config::has_reference_value() == run.exact.has_value(), "reference availability mismatch");
    require(wop::config::has_reference_value() == run.abs_error.has_value(), "abs_error availability mismatch");
}

void test_run_config_wos_matches_direct_solver_statistics_without_reference() {
    const auto path = write_temp_config("runner_wos", cube_config_yaml("wos"));
    const auto cfg = wop::config::load_config_file(path);
    const auto run = wop::config::run_config(cfg);

    std::vector<wop::geometry::Plane> planes;
    planes.reserve(cfg.geometry.planes.size());
    for (const auto& plane : cfg.geometry.planes) {
        planes.push_back(wop::geometry::Plane{plane.p, plane.nu});
    }
    auto poly = wop::geometry::build_polyhedron_from_planes(planes);
    poly = wop::geometry::orient_normals(poly, cfg.geometry.interior_point);

    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return wop::config::boundary_value(y);
    };

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
    require(wop::config::has_reference_value() == run.exact.has_value(), "reference availability mismatch for WoS");
    require(wop::config::has_reference_value() == run.abs_error.has_value(), "abs_error availability mismatch for WoS");
}

void test_cli_accepts_config_and_emits_json(const char* argv0) {
    const auto cfg_path = write_temp_config("cli_wop", cube_config_yaml("wop"));
    const auto cli_path = find_cli_binary(argv0);
    const std::string cmd =
        "cmd /C \"\"" + cli_path.string() + "\" --config \"" + cfg_path.string() + "\" --json\"";
    const std::string out = read_command_output(cmd);

    require(out.find("\"method\":\"wop\"") != std::string::npos, "CLI JSON should contain method");
    require(out.find("\"exact\":") != std::string::npos, "CLI JSON should contain exact");
    require(out.find("\"abs_error\":") != std::string::npos, "CLI JSON should contain abs_error");
}

void test_general_cube_matches_axis_aligned_box_under_identical_wop_parameters() {
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return wop::config::boundary_value(y);
    };

    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    auto axis_box = make_axis_box_polyhedron();
    auto general_cube = make_general_cube_polyhedron();

    wop::rng::Rng rng_box(12345);
    const auto box_result = wop::solver::estimate_wop(
        axis_box,
        x0,
        boundary_f,
        256,
        rng_box,
        std::nullopt,
        std::nullopt,
        1e-14,
        200000,
        0.0,
        1e6,
        wop::solver::RMaxMode::Project,
        3.0);

    wop::rng::Rng rng_general(12345);
    const auto general_result = wop::solver::estimate_wop(
        general_cube,
        x0,
        boundary_f,
        256,
        rng_general,
        std::nullopt,
        std::nullopt,
        1e-14,
        200000,
        0.0,
        1e6,
        wop::solver::RMaxMode::Project,
        3.0);

    require(close(box_result.J, general_result.J, 1e-12), "box/general J mismatch");
    require(close(box_result.S2, general_result.S2, 1e-12), "box/general S2 mismatch");
    require(close(box_result.eps, general_result.eps, 1e-12), "box/general eps mismatch");
    require(close(box_result.mean_steps, general_result.mean_steps, 1e-12), "box/general mean_steps mismatch");
    require(box_result.n_truncated == general_result.n_truncated, "box/general n_truncated mismatch");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        test_problem_functions_are_finite_on_reference_point();
        test_yaml_config_loader_accepts_config_without_boundary_section();
        test_yaml_config_loader_rejects_legacy_boundary_section();
        test_yaml_config_loader_allows_optional_command_metadata();
        test_yaml_config_loader_rejects_missing_interior_point();
        test_example_configs_are_tuned_for_fast_runs(argc > 0 ? argv[0] : "");
        test_run_config_wop_matches_direct_solver_statistics();
        test_run_config_wos_matches_direct_solver_statistics_without_reference();
        test_cli_accepts_config_and_emits_json(argc > 0 ? argv[0] : "");
        test_general_cube_matches_axis_aligned_box_under_identical_wop_parameters();
        std::cout << "wop_config_runtime_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_config_runtime_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
