#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "wop/math/vec3.hpp"
#include "wop/sampling/sampling.hpp"

namespace {

void require(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

wop::math::Vec3 normalize(const wop::math::Vec3& v) {
    return wop::math::normalized(v);
}

double quantile_linear(std::vector<double> values, double p) {
    if (values.empty()) {
        throw std::invalid_argument("quantile_linear requires non-empty values.");
    }
    if (p <= 0.0) {
        return *std::min_element(values.begin(), values.end());
    }
    if (p >= 1.0) {
        return *std::max_element(values.begin(), values.end());
    }

    std::sort(values.begin(), values.end());
    const double pos = p * static_cast<double>(values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);
    return values[lo] * (1.0 - frac) + values[hi] * frac;
}

void test_tangent_sampler_unit_and_orthogonality() {
    wop::rng::Rng rng(20260225);
    const std::vector<wop::math::Vec3> normals = {
        wop::math::Vec3{0.0, 0.0, 1.0},
        normalize(wop::math::Vec3{1.0, -2.0, 0.5}),
        normalize(wop::math::Vec3{-0.7, 0.1, 2.0}),
    };

    for (const auto& nu : normals) {
        for (int i = 0; i < 2000; ++i) {
            const auto w = wop::sampling::sample_tangent_direction(nu, rng);
            require(std::abs(wop::math::norm(w) - 1.0) <= 1e-12, "tangent sample norm mismatch");
            require(std::abs(wop::math::dot(w, nu)) <= 1e-12, "tangent sample orthogonality mismatch");
        }
    }
}

void test_tangent_sampler_second_moment() {
    wop::rng::Rng rng(20260225);
    const auto nu = normalize(wop::math::Vec3{1.0, -2.0, 0.5});

    constexpr int n = 50000;
    std::array<double, 3> mean{0.0, 0.0, 0.0};
    std::array<std::array<double, 3>, 3> second{};

    for (int k = 0; k < n; ++k) {
        const auto w = wop::sampling::sample_tangent_direction(nu, rng);
        const std::array<double, 3> c = {w.x, w.y, w.z};
        for (int i = 0; i < 3; ++i) {
            mean[i] += c[i];
            for (int j = 0; j < 3; ++j) {
                second[i][j] += c[i] * c[j];
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        mean[i] /= static_cast<double>(n);
        require(std::abs(mean[i]) <= 1.5e-2, "tangent sample mean mismatch");
        for (int j = 0; j < 3; ++j) {
            second[i][j] /= static_cast<double>(n);
        }
    }

    const std::array<double, 3> nu_c = {nu.x, nu.y, nu.z};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const double target = 0.5 * ((i == j ? 1.0 : 0.0) - nu_c[i] * nu_c[j]);
            require(std::abs(second[i][j] - target) <= 2.5e-2, "tangent second moment mismatch");
        }
    }
}

double theoretical_quantile_r(double p, double h) {
    return h * std::sqrt(1.0 / ((1.0 - p) * (1.0 - p)) - 1.0);
}

void test_plane_sampler_quantiles() {
    wop::rng::Rng rng(20260216);
    const double h = 1.7;
    const wop::math::Vec3 x{0.0, 0.0, h};
    const wop::math::Vec3 nu{0.0, 0.0, 1.0};
    const double b = 0.0;

    constexpr int n = 60000;
    std::vector<double> radial;
    radial.reserve(n);
    for (int k = 0; k < n; ++k) {
        const auto hit = wop::sampling::sample_hit_on_plane_from_point(x, nu, b, rng);
        radial.push_back(std::hypot(hit.y.x, hit.y.y));
    }

    const std::array<double, 4> probs = {0.25, 0.5, 0.75, 0.9};
    for (double p : probs) {
        const double empirical = quantile_linear(radial, p);
        const double theoretical = theoretical_quantile_r(p, h);
        const double tol = 0.06 * theoretical + 0.02;
        require(std::abs(empirical - theoretical) <= tol, "plane sampler quantile mismatch");
    }
}

}  // namespace

int main() {
    try {
        test_tangent_sampler_unit_and_orthogonality();
        test_tangent_sampler_second_moment();
        test_plane_sampler_quantiles();
        std::cout << "wop_sampling_stats_tests: OK\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "wop_sampling_stats_tests: FAIL: " << ex.what() << "\n";
        return 1;
    }
}
