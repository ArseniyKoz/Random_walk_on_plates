#pragma once

#include <optional>
#include <string>

#include "wop/config/config.hpp"
#include "wop/estimation/estimation.hpp"

namespace wop::config {

struct ConfigRunResult {
    Method method = Method::Wop;
    math::Vec3 x0{};
    estimation::EstimateResult estimate;
    std::optional<double> exact = std::nullopt;
    std::optional<double> abs_error = std::nullopt;
};

ConfigRunResult run_config(const RuntimeConfig& config);
std::string format_json_result(const ConfigRunResult& result);
std::string format_text_result(const ConfigRunResult& result);
const char* method_to_string(Method method) noexcept;

}  // namespace wop::config
