#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>

#include "wop/math/vec3.hpp"
#include "wop/rng/rng.hpp"
#include "wop/solver/wos_box_solver.hpp"

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

}  // namespace

int main() {
    try {
        test_trajectory_stops_on_boundary_delta_neighborhood();
        test_weighted_estimator_is_not_constant_one_for_far_start();
        test_harmonic_reference_solution_is_matched();
        std::cout << "wos_lecture_external_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wos_lecture_external_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
