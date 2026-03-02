#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

#include "wop/estimation/estimation.hpp"
#include "wop/geometry/box.hpp"
#include "wop/solver/wos_box_solver.hpp"

namespace {

void require(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

bool close(double a, double b, double atol = 1e-12) {
    return std::abs(a - b) <= atol;
}

void test_geometry() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    auto poly = wop::geometry::make_axis_aligned_box(mn, mx);
    poly = wop::geometry::orient_normals(poly, wop::math::Vec3{0.0, 0.0, 0.0});

    require(poly.num_faces() == 6, "box must have six faces");
    require(poly.is_inside_or_on(wop::math::Vec3{0.0, 0.0, 0.0}), "origin must be inside");
    require(!poly.is_inside_or_on(wop::math::Vec3{3.0, 0.0, 0.0}), "outside point must be outside");

    const double d = wop::geometry::distance_to_box(wop::math::Vec3{3.0, 0.0, 0.0}, mn, mx);
    require(close(d, 2.0), "distance_to_box mismatch");

    const auto y1 = wop::geometry::closest_point_on_box_boundary(wop::math::Vec3{3.0, 0.2, -0.4}, mn, mx);
    require(close(y1.x, 1.0) && close(y1.y, 0.2) && close(y1.z, -0.4), "closest boundary for outside point mismatch");

    const auto y2 = wop::geometry::closest_point_on_box_boundary(wop::math::Vec3{0.2, 0.3, 0.4}, mn, mx);
    const bool on_boundary = close(std::abs(y2.x), 1.0) || close(std::abs(y2.y), 1.0) || close(std::abs(y2.z), 1.0);
    require(on_boundary, "closest boundary for inside point must lie on boundary");
}

void test_estimation_aggregator() {
    std::vector<wop::estimation::TrajectoryResult> tr = {
        {1.0, 2, "hit_face"},
        {2.0, 3, "timeout"},
        {3.0, 4, "escaped"},
    };
    std::size_t idx = 0;
    const auto result = wop::estimation::estimate_from_trajectories(3, [&]() {
        return tr[idx++];
    });

    require(close(result.J, 2.0), "J mismatch");
    require(close(result.S2, 2.0 / 3.0), "S2 mismatch");
    require(close(result.eps, std::sqrt(2.0)), "eps mismatch");
    require(result.n_total == 3, "n_total mismatch");
    require(result.n_truncated == 2, "n_truncated mismatch");
    require(close(result.mean_steps, 3.0), "mean_steps mismatch");
}

void test_wos_smoke() {
    const wop::math::Vec3 mn{-1.0, -1.0, -1.0};
    const wop::math::Vec3 mx{1.0, 1.0, 1.0};
    const wop::math::Vec3 x0{3.0, 0.0, 0.0};
    const wop::math::Vec3 a{0.2, -0.1, 0.3};

    auto boundary_f = [&](const wop::math::Vec3& y, std::optional<int>) {
        return 1.0 / wop::math::norm(y - a);
    };

    wop::rng::Rng rng(12345);
    const auto res = wop::solver::estimate_wos_box(x0, mn, mx, boundary_f, 512, rng, 1e-3, 200000, 0.0, 1e6);

    require(res.n_total == 512, "wos n_total mismatch");
    require(res.n_truncated >= 0 && res.n_truncated <= res.n_total, "wos n_truncated out of range");
    require(res.S2 >= 0.0, "wos S2 must be non-negative");
    require(close(res.eps, 3.0 * std::sqrt(res.S2 / static_cast<double>(res.n_total)), 1e-12), "wos eps mismatch");
}

}  // namespace

int main() {
    try {
        test_geometry();
        test_estimation_aggregator();
        test_wos_smoke();
        std::cout << "wop_core_smoke_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_core_smoke_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
