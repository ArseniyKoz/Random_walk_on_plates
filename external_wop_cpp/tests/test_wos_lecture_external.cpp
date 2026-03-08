#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>

#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"
#include "wop/solver/wos_box_solver.hpp"
#include "wop/solver/wos_solver.hpp"
#include "wop/geometry/box.hpp"

namespace {

void require(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

void test_trajectory_stops_on_boundary_delta_neighborhood() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    const auto boundary_f = [](const wop::math::Vec3&, std::optional<int>) { return 1.0; };
    wop::rng::Rng rng(20260304);

    const auto tr = wop::solver::trace_wos_box_trajectory(
        x0,
        mn,
        mx,
        boundary_f,
        rng,
        1e-3,
        2.0,
        4.0,
        100000,
        0.0);

    require(tr.status == "hit_face", "lecture WoS should stop via boundary neighborhood hit");
    require(tr.steps > 0, "trajectory should perform at least one step");
}

void test_near_boundary_returns_boundary_value() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{1.0 + 5e-4, 0.2, -0.3};
    const auto boundary_f = [](const wop::math::Vec3& y, std::optional<int>) { return y.x; };
    wop::rng::Rng rng(20260307);

    const auto tr = wop::solver::trace_wos_box_trajectory(
        x0,
        mn,
        mx,
        boundary_f,
        rng,
        1e-3,
        1.0,
        2.0,
        1000,
        0.0);

    require(tr.status == "hit_face", "near-boundary WoS must terminate as hit_face");
    require(tr.steps == 1, "near-boundary WoS must terminate on first step");
    require(std::abs(tr.value - 1.0) <= 1e-12, "boundary value must be evaluated at projected boundary point");
}

void test_weighted_estimator_is_not_constant_one_for_far_start() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{10.0, 0.0, 0.0};
    const auto boundary_f = [](const wop::math::Vec3&, std::optional<int>) { return 1.0; };
    wop::rng::Rng rng(20260305);

    const auto res = wop::solver::estimate_wos_box(
        x0,
        mn,
        mx,
        boundary_f,
        2000,
        rng,
        1e-3,
        2.0,
        4.0,
        200000,
        0.0);

    require(res.n_truncated >= 0 && res.n_truncated <= res.n_total, "n_truncated out of range");
    require(res.J < 0.95, "for far start and constant boundary, weighted lecture estimator must be below 1");
}

void test_harmonic_reference_solution_is_matched() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.2, -0.1, 0.3};
    const auto exact_u = [&](const wop::math::Vec3& x) { return 1.0 / wop::math::norm(x - a); };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };
    wop::rng::Rng rng(20260306);

    const auto res = wop::solver::estimate_wos_box(
        x0,
        mn,
        mx,
        boundary_f,
        3000,
        rng,
        1e-3,
        2.0,
        4.0,
        200000,
        0.0);

    const double exact = exact_u(x0);
    const double err = std::abs(res.J - exact);
    require(err <= 2.5 * res.eps + 1e-2, "lecture WoS should match harmonic reference within MC tolerance");
}

void test_harmonic_reference_solution_is_matched_for_rectangular_box() {
    const wop::math::Vec3 mn{-2.0, -1.0, -0.5};
    const wop::math::Vec3 mx{1.0, 2.5, 1.5};
    const wop::math::Vec3 x0{4.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.0, 0.3, 0.2};
    const auto exact_u = [&](const wop::math::Vec3& x) { return 1.0 / wop::math::norm(x - a); };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return exact_u(y); };
    wop::rng::Rng rng(20260308);

    const auto res = wop::solver::estimate_wos_box(
        x0,
        mn,
        mx,
        boundary_f,
        3000,
        rng,
        1e-3,
        1.0,
        2.0,
        200000,
        0.0);

    const double exact = exact_u(x0);
    const double err = std::abs(res.J - exact);
    require(err <= 2.5 * res.eps + 1.5e-2, "lecture WoS should match harmonic reference for rectangular box");
}

void test_general_wos_matches_box_wrapper_for_unit_box() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.2, -0.1, 0.3};
    auto poly = wop::geometry::make_axis_aligned_box(mn, mx);
    poly = wop::geometry::orient_normals(poly, wop::math::Vec3{0.0, 0.0, 0.0});
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) { return 1.0 / wop::math::norm(y - a); };

    wop::rng::Rng rng_poly(20260309);
    const auto res_poly = wop::solver::estimate_wos(
        poly,
        x0,
        boundary_f,
        3000,
        rng_poly,
        1e-3,
        1.0,
        2.0,
        200000,
        0.0);

    wop::rng::Rng rng_box(20260309);
    const auto res_box = wop::solver::estimate_wos_box(
        x0,
        mn,
        mx,
        boundary_f,
        3000,
        rng_box,
        1e-3,
        1.0,
        2.0,
        200000,
        0.0);

    const double diff = std::abs(res_poly.J - res_box.J);
    const double se_poly = std::sqrt(res_poly.S2 / static_cast<double>(res_poly.n_total));
    const double se_box = std::sqrt(res_box.S2 / static_cast<double>(res_box.n_total));
    const double combined = std::sqrt(se_poly * se_poly + se_box * se_box);
    require(diff <= 4.0 * combined + 2e-3, "general WoS and box wrapper must agree statistically");
}

}  // namespace

int main() {
    try {
        test_trajectory_stops_on_boundary_delta_neighborhood();
        test_near_boundary_returns_boundary_value();
        test_weighted_estimator_is_not_constant_one_for_far_start();
        test_harmonic_reference_solution_is_matched();
        test_harmonic_reference_solution_is_matched_for_rectangular_box();
        test_general_wos_matches_box_wrapper_for_unit_box();
        std::cout << "wos_lecture_external_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wos_lecture_external_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
