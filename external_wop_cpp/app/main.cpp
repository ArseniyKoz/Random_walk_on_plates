#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wop/config/config.hpp"
#include "wop/config/run_config.hpp"
#include "wop/geometry/box.hpp"
#include "wop/solver/wop_solver.hpp"
#include "wop/solver/wos_solver.hpp"

namespace {

struct CliArgs {
    std::optional<std::string> config_path = std::nullopt;
    std::string example = "box";
    std::string method = "wop";
    wop::math::Vec3 x0{3.0, 0.0, 0.0};
    int n = 50000;
    std::uint64_t seed = 12345;
    int max_steps = 1'000'000;
    double delta = 1e-3;
    double rho_scale = 1.0;
    double rho1_scale = 2.0;
    std::optional<double> r_max = 1e6;
    wop::solver::RMaxMode r_max_mode = wop::solver::RMaxMode::Project;
    double r_max_factor = 3.0;
    bool json = false;
    bool legacy_mode_args_used = false;
};

[[noreturn]] void throw_usage_error(const std::string& msg) {
    throw std::invalid_argument(msg);
}

wop::math::Vec3 parse_vec3_arg(const std::string& text) {
    std::istringstream iss(text);
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    std::string extra;
    if (!(iss >> x >> y >> z) || (iss >> extra)) {
        throw std::invalid_argument("Invalid --x0 format, expected: \"3 0 0\".");
    }
    return wop::math::Vec3{x, y, z};
}

wop::solver::RMaxMode parse_r_max_mode(const std::string& text) {
    if (text == "escape") {
        return wop::solver::RMaxMode::Escape;
    }
    if (text == "project") {
        return wop::solver::RMaxMode::Project;
    }
    throw std::invalid_argument("Invalid --r-max-mode, expected one of: escape, project.");
}

CliArgs parse_args(int argc, char** argv) {
    CliArgs args;

    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto require_value = [&](const std::string& opt) -> std::string {
            if (i + 1 >= argc) {
                throw_usage_error("Missing value for " + opt + ".");
            }
            return argv[++i];
        };

        if (key == "--example") {
            args.legacy_mode_args_used = true;
            args.example = require_value("--example");
        } else if (key == "--method") {
            args.legacy_mode_args_used = true;
            args.method = require_value("--method");
        } else if (key == "--x0") {
            args.legacy_mode_args_used = true;
            args.x0 = parse_vec3_arg(require_value("--x0"));
        } else if (key == "--n") {
            args.legacy_mode_args_used = true;
            args.n = std::stoi(require_value("--n"));
        } else if (key == "--seed") {
            args.legacy_mode_args_used = true;
            args.seed = static_cast<std::uint64_t>(std::stoull(require_value("--seed")));
        } else if (key == "--max-steps") {
            args.legacy_mode_args_used = true;
            args.max_steps = std::stoi(require_value("--max-steps"));
        } else if (key == "--delta") {
            args.legacy_mode_args_used = true;
            args.delta = std::stod(require_value("--delta"));
        } else if (key == "--rho-scale") {
            args.legacy_mode_args_used = true;
            args.rho_scale = std::stod(require_value("--rho-scale"));
        } else if (key == "--rho1-scale") {
            args.legacy_mode_args_used = true;
            args.rho1_scale = std::stod(require_value("--rho1-scale"));
        } else if (key == "--r-max") {
            args.legacy_mode_args_used = true;
            const double value = std::stod(require_value("--r-max"));
            args.r_max = (value > 0.0) ? std::optional<double>(value) : std::nullopt;
        } else if (key == "--r-max-mode") {
            args.legacy_mode_args_used = true;
            args.r_max_mode = parse_r_max_mode(require_value("--r-max-mode"));
        } else if (key == "--r-max-factor") {
            args.legacy_mode_args_used = true;
            args.r_max_factor = std::stod(require_value("--r-max-factor"));
        } else if (key == "--config") {
            args.config_path = require_value("--config");
        } else if (key == "--json") {
            args.json = true;
        } else if (key == "--help" || key == "-h") {
            throw_usage_error("");
        } else {
            throw_usage_error("Unknown argument: " + key);
        }
    }

    if (args.config_path.has_value()) {
        if (args.legacy_mode_args_used) {
            throw_usage_error("--config cannot be combined with legacy example arguments.");
        }
        return args;
    }

    if (args.n <= 0) {
        throw_usage_error("--n must be positive.");
    }
    if (args.max_steps <= 0) {
        throw_usage_error("--max-steps must be positive.");
    }
    if (args.method != "wop" && args.method != "wos") {
        throw_usage_error("--method must be one of: wop, wos.");
    }
    if (args.delta <= 0.0) {
        throw_usage_error("--delta must be positive.");
    }
    if (args.rho_scale <= 0.0) {
        throw_usage_error("--rho-scale must be positive.");
    }
    if (args.rho1_scale <= 1.0) {
        throw_usage_error("--rho1-scale must be greater than 1.0.");
    }
    if (args.r_max_factor <= 1.0) {
        throw_usage_error("--r-max-factor must be greater than 1.0.");
    }
    if (args.example != "box") {
        throw_usage_error("Unsupported example: " + args.example);
    }

    return args;
}

void print_usage() {
    std::cout << "Usage: wop_cli --config config.yaml [--json]\n"
              << "       wop_cli [--method wop|wos] [--example box] [--x0 \"3 0 0\"] [--n 50000] [--seed 12345]\n"
              << "               [--max-steps 1000000] [--json]\n"
              << "               WOP args: [--r-max 1e6] [--r-max-mode escape|project] [--r-max-factor 3.0]\n"
              << "               WoS args: [--delta 1e-3] [--rho-scale 1.0] [--rho1-scale 2.0]\n";
}

int run_config_mode(const CliArgs& args) {
    const wop::config::RuntimeConfig config = wop::config::load_config_file(*args.config_path);
    const wop::config::ConfigRunResult result = wop::config::run_config(config);
    if (args.json) {
        std::cout << wop::config::format_json_result(result);
    } else {
        std::cout << wop::config::format_text_result(result);
    }
    return 0;
}

int run_box_example(const CliArgs& args) {
    const wop::math::Vec3 box_min{-1.0, -1.0, -1.0};
    const wop::math::Vec3 box_max{1.0, 1.0, 1.0};
    const wop::math::Vec3 interior{0.0, 0.0, 0.0};

    auto poly = wop::geometry::make_axis_aligned_box(box_min, box_max);
    poly = wop::geometry::orient_normals(poly, interior);

    const wop::math::Vec3 a{0.2, -0.1, 0.3};
    const auto exact_u = [&](const wop::math::Vec3& x) {
        return 1.0 / wop::math::norm(x - a);
    };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return exact_u(y);
    };

    wop::rng::Rng rng(args.seed);
    wop::estimation::EstimateResult result{};
    if (args.method == "wop") {
        result = wop::solver::estimate_wop(
            poly,
            args.x0,
            boundary_f,
            args.n,
            rng,
            std::nullopt,
            std::nullopt,
            1e-14,
            args.max_steps,
            0.0,
            args.r_max,
            args.r_max_mode,
            args.r_max_factor);
    } else {
        result = wop::solver::estimate_wos(
            poly,
            args.x0,
            boundary_f,
            args.n,
            rng,
            args.delta,
            args.rho_scale,
            args.rho1_scale,
            args.max_steps,
            0.0);
    }

    const double exact = exact_u(args.x0);
    const double abs_err = std::abs(result.J - exact);

    if (args.json) {
        std::cout << std::setprecision(17)
                  << "{"
                  << "\"method\":\"" << args.method << "\","
                  << "\"x0\":[" << args.x0.x << "," << args.x0.y << "," << args.x0.z << "],"
                  << "\"n_total\":" << result.n_total << ","
                  << "\"n_truncated\":" << result.n_truncated << ","
                  << "\"J\":" << result.J << ","
                  << "\"exact\":" << exact << ","
                  << "\"abs_error\":" << abs_err << ","
                  << "\"S2\":" << result.S2 << ","
                  << "\"eps\":" << result.eps << ","
                  << "\"mean_steps\":" << result.mean_steps
                  << "}\n";
        return 0;
    }

    std::cout << "method: " << args.method << "\n";
    std::cout << "x0: [" << args.x0.x << ", " << args.x0.y << ", " << args.x0.z << "]\n";
    std::cout << "n_total: " << result.n_total << "\n";
    std::cout << "n_truncated: " << result.n_truncated << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "J: " << result.J << "\n";
    std::cout << "exact: " << exact << "\n";
    std::cout << "abs_error: " << abs_err << "\n";
    std::cout << "S2: " << result.S2 << "\n";
    std::cout << "eps: " << result.eps << "\n";
    std::cout << std::setprecision(2) << "mean_steps: " << result.mean_steps << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliArgs args = parse_args(argc, argv);
        if (args.config_path.has_value()) {
            return run_config_mode(args);
        }
        return run_box_example(args);
    } catch (const std::invalid_argument& ex) {
        if (std::string(ex.what()).empty()) {
            print_usage();
            return 0;
        }
        std::cerr << "Error: " << ex.what() << "\n";
        print_usage();
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << "Fatal: " << ex.what() << "\n";
        return 1;
    }
}
