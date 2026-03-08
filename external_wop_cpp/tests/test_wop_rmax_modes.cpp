#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>

#include "wop/geometry/box.hpp"
#include "wop/math/vec3.hpp"
#include "wop/solver/wop_solver.hpp"
#include "wop/solver/wop_solver_internal.hpp"

namespace {

void require(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

bool close(double a, double b, double atol = 1e-12) {
    return std::abs(a - b) <= atol;
}

wop::geometry::Polyhedron make_unit_box_poly() {
    auto poly = wop::geometry::make_axis_aligned_box(
        wop::math::Vec3{-1.0, -1.0, -1.0},
        wop::math::Vec3{1.0, 1.0, 1.0});
    return wop::geometry::orient_normals(poly, wop::math::Vec3{0.0, 0.0, 0.0});
}

wop::solver::BoundaryFunc make_constant_boundary(double value) {
    return [value](const wop::math::Vec3&, std::optional<int>) { return value; };
}

void test_escape_mode_keeps_legacy_escape_behavior() {
    const auto poly = make_unit_box_poly();
    const auto boundary_f = make_constant_boundary(1.0);
    wop::rng::Rng rng(20260303);

    const auto tr = wop::solver::trace_wop_trajectory(
        poly,
        wop::math::Vec3{3.0, 0.0, 0.0},
        boundary_f,
        rng,
        1e-12,
        1e-12,
        1e-3,
        1000,
        0.0,
        2.0,
        wop::solver::RMaxMode::Escape,
        3.0);

    require(tr.status == "escaped", "escape mode must keep escaped truncation");
    require(tr.steps == 0, "escape mode with |x0|>=r_max must truncate at step 0");
}

void test_project_mode_does_not_escape_on_rmax_crossing() {
    const auto poly = make_unit_box_poly();
    const auto boundary_f = make_constant_boundary(1.0);
    wop::rng::Rng rng(20260303);

    const auto tr = wop::solver::trace_wop_trajectory(
        poly,
        wop::math::Vec3{3.0, 0.0, 0.0},
        boundary_f,
        rng,
        1e-12,
        1e-12,
        1e-3,
        1000,
        0.0,
        2.0,
        wop::solver::RMaxMode::Project,
        3.0);

    require(tr.status != "escaped", "project mode must not terminate as escaped on r_max crossing");
}

void test_project_mode_autormax_is_enabled_when_rmax_is_none() {
    const auto poly = make_unit_box_poly();
    const auto boundary_f = make_constant_boundary(1.0);
    wop::rng::Rng rng(20260303);

    const auto tr = wop::solver::trace_wop_trajectory(
        poly,
        wop::math::Vec3{4.0, 0.5, 0.0},
        boundary_f,
        rng,
        1e-12,
        1e-12,
        1e-3,
        1000,
        0.0,
        std::nullopt,
        wop::solver::RMaxMode::Project,
        3.0);

    require(tr.status != "escaped", "project mode with auto-r_max should avoid escaped on sphere crossing");
}

void test_manual_rmax_overrides_auto_rmax() {
    const auto poly = make_unit_box_poly();
    const wop::math::Vec3 x0{4.0, 0.0, 0.0};

    const auto cfg_auto = wop::solver::detail::resolve_r_max_projection(
        poly,
        x0,
        std::nullopt,
        wop::solver::RMaxMode::Project,
        3.0);
    const auto cfg_manual = wop::solver::detail::resolve_r_max_projection(
        poly,
        x0,
        2.5,
        wop::solver::RMaxMode::Project,
        3.0);

    require(cfg_auto.enabled, "auto project config must be enabled");
    require(cfg_manual.enabled, "manual project config must be enabled");
    require(close(cfg_manual.rho, 2.5), "manual r_max must override manual rho");
    require(close(cfg_manual.rho1, 7.5), "manual rho1 must be scaled by r_max_factor");
    require(cfg_auto.rho > cfg_manual.rho, "auto and manual rho should differ in this configuration");
}

void test_project_mode_is_not_systematically_biased_for_exact_harmonic_solution() {
    const auto poly = make_unit_box_poly();
    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.2, -0.1, 0.3};
    const auto exact_u = [&](const wop::math::Vec3& x) { return 1.0 / wop::math::norm(x - a); };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };

    wop::rng::Rng rng_escape(20260304);
    const auto res_escape = wop::solver::estimate_wop(
        poly,
        x0,
        boundary_f,
        10000,
        rng_escape,
        std::nullopt,
        std::nullopt,
        1e-14,
        400000,
        0.0,
        1e6,
        wop::solver::RMaxMode::Escape,
        3.0);

    wop::rng::Rng rng_project(20260304);
    const auto res_project = wop::solver::estimate_wop(
        poly,
        x0,
        boundary_f,
        10000,
        rng_project,
        std::nullopt,
        std::nullopt,
        1e-14,
        400000,
        0.0,
        2.5,
        wop::solver::RMaxMode::Project,
        2.0);

    const double exact = exact_u(x0);
    const double err_escape = std::abs(res_escape.J - exact);
    const double err_project = std::abs(res_project.J - exact);

    require(err_project < 0.15, "project mode should not have large systematic bias");
    require(std::abs(res_project.J - res_escape.J) < 0.15, "project and escape estimates should be comparable");
    require(res_project.n_total == 10000, "project n_total mismatch");
    require(res_escape.n_total == 10000, "escape n_total mismatch");
    require(err_escape < 0.1, "escape mode baseline unexpectedly inaccurate");
}

void test_project_mode_handles_nonzero_uinf_shift() {
    const auto poly = make_unit_box_poly();
    const wop::math::Vec3 x0{6.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.2, -0.1, 0.3};
    constexpr double u_inf = 2.0;
    const auto exact_u = [&](const wop::math::Vec3& x) { return u_inf + 1.0 / wop::math::norm(x - a); };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };

    wop::rng::Rng rng_project(20260305);
    const auto res_project = wop::solver::estimate_wop(
        poly,
        x0,
        boundary_f,
        12000,
        rng_project,
        std::nullopt,
        std::nullopt,
        1e-14,
        500000,
        u_inf,
        2.5,
        wop::solver::RMaxMode::Project,
        2.0);

    const double exact = exact_u(x0);
    const double err = std::abs(res_project.J - exact);

    require(err < 0.2, "project mode with non-zero u_inf should stay close to exact value");
    require(res_project.n_total == 12000, "project n_total mismatch for non-zero u_inf case");
}

void test_invalid_rmax_factor_raises() {
    const auto poly = make_unit_box_poly();
    const auto boundary_f = make_constant_boundary(1.0);
    wop::rng::Rng rng(20260303);

    bool caught = false;
    try {
        static_cast<void>(wop::solver::trace_wop_trajectory(
            poly,
            wop::math::Vec3{3.0, 0.0, 0.0},
            boundary_f,
            rng,
            1e-12,
            1e-12,
            1e-3,
            1000,
            0.0,
            std::nullopt,
            wop::solver::RMaxMode::Project,
            1.0));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "r_max_factor <= 1 must be rejected");
}

}  // namespace

int main() {
    try {
        test_escape_mode_keeps_legacy_escape_behavior();
        test_project_mode_does_not_escape_on_rmax_crossing();
        test_project_mode_autormax_is_enabled_when_rmax_is_none();
        test_manual_rmax_overrides_auto_rmax();
        test_project_mode_is_not_systematically_biased_for_exact_harmonic_solution();
        test_project_mode_handles_nonzero_uinf_shift();
        test_invalid_rmax_factor_raises();
        std::cout << "wop_rmax_modes_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_rmax_modes_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
