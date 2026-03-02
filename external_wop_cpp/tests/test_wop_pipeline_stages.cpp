#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

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

void test_scan_distances_outside_and_hit_flags() {
    const std::vector<double> values{-0.2, 0.3, -0.1, 0.05};
    const auto scan = wop::solver::detail::scan_distances(values, 0.0, 3, 1e-3);

    require(scan.any_outside, "scan.any_outside mismatch");
    require(!scan.all_inside, "scan.all_inside mismatch");
    require(scan.argmin_outside == 3, "scan.argmin_outside mismatch");
    require(scan.argmin_abs == 3, "scan.argmin_abs mismatch");
    require(!scan.hit_target_face, "scan.hit_target_face mismatch");
}

void test_scan_distances_all_inside_and_hit_target() {
    const std::vector<double> values{-0.2, -1e-5, -0.1, -1e-6};
    const auto scan = wop::solver::detail::scan_distances(values, 0.0, 1, 1e-4);

    require(!scan.any_outside, "scan.any_outside expected false");
    require(scan.all_inside, "scan.all_inside expected true");
    require(scan.argmin_outside == 0, "scan.argmin_outside default mismatch");
    require(scan.argmin_abs == 3, "scan.argmin_abs mismatch");
    require(scan.hit_target_face, "scan.hit_target_face expected true");
}

void test_sample_plane_radius_domain_and_finiteness() {
    wop::rng::Rng rng(20260302);
    for (int i = 0; i < 2000; ++i) {
        const double rho = wop::solver::detail::sample_plane_radius(0.7, rng, 1e-14);
        require(std::isfinite(rho), "sample_plane_radius must be finite");
        require(rho >= 0.0, "sample_plane_radius must be non-negative");
    }

    bool caught = false;
    try {
        static_cast<void>(wop::solver::detail::sample_plane_radius(0.0, rng, 1e-14));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "sample_plane_radius should reject non-positive d");

    caught = false;
    try {
        static_cast<void>(wop::solver::detail::sample_plane_radius(1.0, rng, 0.0));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "sample_plane_radius should reject non-positive min_one_minus_alpha");
}

void test_unit_cube_n1000_matches_exact_solution() {
    const wop::math::Vec3 mn{0.0, 0.0, 0.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 interior{0.5, 0.5, 0.5};

    auto poly = wop::geometry::make_axis_aligned_box(mn, mx);
    poly = wop::geometry::orient_normals(poly, interior);

    const wop::math::Vec3 a{0.25, 0.35, 0.45};
    const wop::math::Vec3 x0{2.0, 0.5, 0.5};

    const auto exact_u = [&](const wop::math::Vec3& x) {
        return 1.0 / wop::math::norm(x - a);
    };
    const auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return exact_u(y);
    };

    wop::rng::Rng rng(20260302);
    const auto result = wop::solver::estimate_wop(
        poly,
        x0,
        boundary_f,
        1000,
        rng,
        std::nullopt,
        std::nullopt,
        1e-14,
        200000,
        0.0,
        1e6);

    require(result.n_total == 1000, "n_total mismatch");
    require(result.n_truncated >= 0 && result.n_truncated <= result.n_total, "n_truncated out of range");
    require(result.S2 >= 0.0, "S2 must be non-negative");
    require(close(result.eps, 3.0 * std::sqrt(result.S2 / static_cast<double>(result.n_total))), "eps mismatch");

    const double exact = exact_u(x0);
    const double err = std::abs(result.J - exact);
    const double tol = 2.0 * result.eps + 1e-2;
    require(err <= tol, "unit cube n=1000 estimator mismatch");
}

}  // namespace

int main() {
    try {
        test_scan_distances_outside_and_hit_flags();
        test_scan_distances_all_inside_and_hit_target();
        test_sample_plane_radius_domain_and_finiteness();
        test_unit_cube_n1000_matches_exact_solution();
        std::cout << "wop_pipeline_stages_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_pipeline_stages_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
