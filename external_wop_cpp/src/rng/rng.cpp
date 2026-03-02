#include "wop/rng/rng.hpp"

namespace wop::rng {

Rng::Rng(std::uint64_t seed)
    : engine_(seed), uniform_dist_(0.0, 1.0), normal_dist_(0.0, 1.0) {}

double Rng::uniform01() {
    return uniform_dist_(engine_);
}

double Rng::normal01() {
    return normal_dist_(engine_);
}

math::Vec3 Rng::normal3() {
    return math::Vec3{normal01(), normal01(), normal01()};
}

}  // namespace wop::rng
