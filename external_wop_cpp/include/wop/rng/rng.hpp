#pragma once

#include <cstdint>
#include <random>

#include "wop/math/vec3.hpp"

namespace wop::rng {

class Rng {
public:
    explicit Rng(std::uint64_t seed);

    double uniform01();
    double normal01();
    math::Vec3 normal3();

private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> uniform_dist_;
    std::normal_distribution<double> normal_dist_;
};

}  // namespace wop::rng
