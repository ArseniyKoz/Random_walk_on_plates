#pragma once

#include <string>
#include <vector>

#include "wop/math/vec3.hpp"

namespace wop::config {

class CompiledExpression {
public:
    enum class OpCode {
        PushConstant,
        PushX,
        PushY,
        PushZ,
        UnaryPlus,
        UnaryMinus,
        Add,
        Subtract,
        Multiply,
        Divide,
        Power,
        Sqrt,
        Exp,
        Log,
        Sin,
        Cos,
        Abs,
    };

    struct Instruction {
        OpCode op = OpCode::PushConstant;
        double value = 0.0;
    };

    CompiledExpression() = default;

    double evaluate(const math::Vec3& point) const;
    const std::string& source() const noexcept { return source_; }

private:
    std::vector<Instruction> program_;
    std::string source_;

    CompiledExpression(std::vector<Instruction> program, std::string source);

    friend CompiledExpression compile_expression(const std::string& text);
};

CompiledExpression compile_expression(const std::string& text);

}  // namespace wop::config
