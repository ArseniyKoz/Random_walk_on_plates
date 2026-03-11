#include "wop/config/expression.hpp"

#include <cmath>
#include <cctype>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace wop::config {

namespace {

using OpCode = CompiledExpression::OpCode;
using Instruction = CompiledExpression::Instruction;

enum class TokenType {
    Number,
    Identifier,
    Operator,
    LeftParen,
    RightParen,
};

struct Token {
    TokenType type = TokenType::Number;
    std::string text;
    double number = 0.0;
};

bool is_identifier_start(char ch) {
    return std::isalpha(static_cast<unsigned char>(ch)) != 0 || ch == '_';
}

bool is_identifier_continue(char ch) {
    return std::isalnum(static_cast<unsigned char>(ch)) != 0 || ch == '_';
}

std::vector<Token> tokenize(const std::string& text) {
    std::vector<Token> tokens;
    for (std::size_t i = 0; i < text.size();) {
        const char ch = text[i];
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            ++i;
            continue;
        }

        if (std::isdigit(static_cast<unsigned char>(ch)) != 0 || ch == '.') {
            std::size_t consumed = 0;
            const double value = std::stod(text.substr(i), &consumed);
            tokens.push_back(Token{TokenType::Number, text.substr(i, consumed), value});
            i += consumed;
            continue;
        }

        if (is_identifier_start(ch)) {
            std::size_t j = i + 1;
            while (j < text.size() && is_identifier_continue(text[j])) {
                ++j;
            }
            tokens.push_back(Token{TokenType::Identifier, text.substr(i, j - i), 0.0});
            i = j;
            continue;
        }

        if (ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^') {
            tokens.push_back(Token{TokenType::Operator, std::string(1, ch), 0.0});
            ++i;
            continue;
        }

        if (ch == '(') {
            tokens.push_back(Token{TokenType::LeftParen, "(", 0.0});
            ++i;
            continue;
        }

        if (ch == ')') {
            tokens.push_back(Token{TokenType::RightParen, ")", 0.0});
            ++i;
            continue;
        }

        throw std::invalid_argument("Unsupported token in expression: " + std::string(1, ch));
    }
    return tokens;
}

bool is_function_name(const std::string& text) {
    return text == "sqrt" || text == "exp" || text == "log" || text == "sin" || text == "cos" || text == "abs";
}

Instruction make_identifier_instruction(const std::string& text) {
    if (text == "x") {
        return Instruction{OpCode::PushX, 0.0};
    }
    if (text == "y") {
        return Instruction{OpCode::PushY, 0.0};
    }
    if (text == "z") {
        return Instruction{OpCode::PushZ, 0.0};
    }
    if (text == "pi") {
        return Instruction{OpCode::PushConstant, std::acos(-1.0)};
    }
    if (text == "e") {
        return Instruction{OpCode::PushConstant, std::exp(1.0)};
    }
    throw std::invalid_argument("Unknown identifier in expression: " + text);
}

Instruction make_function_instruction(const std::string& text) {
    if (text == "sqrt") {
        return Instruction{OpCode::Sqrt, 0.0};
    }
    if (text == "exp") {
        return Instruction{OpCode::Exp, 0.0};
    }
    if (text == "log") {
        return Instruction{OpCode::Log, 0.0};
    }
    if (text == "sin") {
        return Instruction{OpCode::Sin, 0.0};
    }
    if (text == "cos") {
        return Instruction{OpCode::Cos, 0.0};
    }
    if (text == "abs") {
        return Instruction{OpCode::Abs, 0.0};
    }
    throw std::invalid_argument("Unknown function in expression: " + text);
}

struct StackItem {
    enum class Kind {
        Operator,
        Function,
        LeftParen,
    };

    Kind kind = Kind::Operator;
    std::string text;
};

struct OperatorInfo {
    int precedence = 0;
    bool right_associative = false;
    Instruction instruction;
};

OperatorInfo get_operator_info(const std::string& op) {
    if (op == "u+") {
        return OperatorInfo{3, true, Instruction{OpCode::UnaryPlus, 0.0}};
    }
    if (op == "u-") {
        return OperatorInfo{3, true, Instruction{OpCode::UnaryMinus, 0.0}};
    }
    if (op == "^") {
        return OperatorInfo{4, true, Instruction{OpCode::Power, 0.0}};
    }
    if (op == "*") {
        return OperatorInfo{2, false, Instruction{OpCode::Multiply, 0.0}};
    }
    if (op == "/") {
        return OperatorInfo{2, false, Instruction{OpCode::Divide, 0.0}};
    }
    if (op == "+") {
        return OperatorInfo{1, false, Instruction{OpCode::Add, 0.0}};
    }
    if (op == "-") {
        return OperatorInfo{1, false, Instruction{OpCode::Subtract, 0.0}};
    }
    throw std::invalid_argument("Unknown operator in expression: " + op);
}

void append_stack_item(std::vector<Instruction>& output, const StackItem& item) {
    if (item.kind == StackItem::Kind::Operator) {
        output.push_back(get_operator_info(item.text).instruction);
        return;
    }
    if (item.kind == StackItem::Kind::Function) {
        output.push_back(make_function_instruction(item.text));
        return;
    }
    throw std::invalid_argument("Mismatched parentheses in expression.");
}

std::vector<Instruction> compile_to_rpn(const std::vector<Token>& tokens) {
    if (tokens.empty()) {
        throw std::invalid_argument("Expression must not be empty.");
    }

    std::vector<Instruction> output;
    std::vector<StackItem> stack;
    bool expect_operand = true;

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const Token& token = tokens[i];

        if (token.type == TokenType::Number) {
            if (!expect_operand) {
                throw std::invalid_argument("Missing operator before number.");
            }
            output.push_back(Instruction{OpCode::PushConstant, token.number});
            expect_operand = false;
            continue;
        }

        if (token.type == TokenType::Identifier) {
            const bool is_function = (i + 1 < tokens.size() && tokens[i + 1].type == TokenType::LeftParen);
            if (is_function) {
                if (!is_function_name(token.text)) {
                    throw std::invalid_argument("Unknown function in expression: " + token.text);
                }
                stack.push_back(StackItem{StackItem::Kind::Function, token.text});
                expect_operand = true;
            } else {
                if (!expect_operand) {
                    throw std::invalid_argument("Missing operator before identifier.");
                }
                output.push_back(make_identifier_instruction(token.text));
                expect_operand = false;
            }
            continue;
        }

        if (token.type == TokenType::LeftParen) {
            if (!expect_operand) {
                throw std::invalid_argument("Implicit multiplication is not supported.");
            }
            stack.push_back(StackItem{StackItem::Kind::LeftParen, "("});
            expect_operand = true;
            continue;
        }

        if (token.type == TokenType::RightParen) {
            if (expect_operand) {
                throw std::invalid_argument("Unexpected ')' in expression.");
            }
            bool found_left = false;
            while (!stack.empty()) {
                const StackItem item = stack.back();
                stack.pop_back();
                if (item.kind == StackItem::Kind::LeftParen) {
                    found_left = true;
                    break;
                }
                append_stack_item(output, item);
            }
            if (!found_left) {
                throw std::invalid_argument("Mismatched parentheses in expression.");
            }
            if (!stack.empty() && stack.back().kind == StackItem::Kind::Function) {
                append_stack_item(output, stack.back());
                stack.pop_back();
            }
            expect_operand = false;
            continue;
        }

        std::string op = token.text;
        if (expect_operand) {
            if (op == "+" || op == "-") {
                op = (op == "+") ? "u+" : "u-";
            } else {
                throw std::invalid_argument("Unexpected binary operator in expression.");
            }
        }

        const OperatorInfo current = get_operator_info(op);
        while (!stack.empty() && stack.back().kind == StackItem::Kind::Operator) {
            const OperatorInfo top = get_operator_info(stack.back().text);
            const bool should_pop =
                (!current.right_associative && current.precedence <= top.precedence) ||
                (current.right_associative && current.precedence < top.precedence);
            if (!should_pop) {
                break;
            }
            append_stack_item(output, stack.back());
            stack.pop_back();
        }

        stack.push_back(StackItem{StackItem::Kind::Operator, op});
        expect_operand = true;
    }

    if (expect_operand) {
        throw std::invalid_argument("Expression ended unexpectedly.");
    }

    while (!stack.empty()) {
        if (stack.back().kind == StackItem::Kind::LeftParen) {
            throw std::invalid_argument("Mismatched parentheses in expression.");
        }
        append_stack_item(output, stack.back());
        stack.pop_back();
    }

    return output;
}

double pop_one(std::vector<double>& stack, std::string_view name) {
    if (stack.empty()) {
        throw std::invalid_argument("Expression stack underflow while evaluating " + std::string(name) + ".");
    }
    const double value = stack.back();
    stack.pop_back();
    return value;
}

void ensure_finite(double value, std::string_view source) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Expression evaluated to a non-finite value: " + std::string(source));
    }
}

}  // namespace

CompiledExpression::CompiledExpression(std::vector<Instruction> program, std::string source)
    : program_(std::move(program)), source_(std::move(source)) {}

double CompiledExpression::evaluate(const math::Vec3& point) const {
    std::vector<double> stack;
    stack.reserve(program_.size());

    for (const Instruction& ins : program_) {
        switch (ins.op) {
            case OpCode::PushConstant:
                stack.push_back(ins.value);
                break;
            case OpCode::PushX:
                stack.push_back(point.x);
                break;
            case OpCode::PushY:
                stack.push_back(point.y);
                break;
            case OpCode::PushZ:
                stack.push_back(point.z);
                break;
            case OpCode::UnaryPlus:
                stack.push_back(pop_one(stack, "unary +"));
                break;
            case OpCode::UnaryMinus:
                stack.push_back(-pop_one(stack, "unary -"));
                break;
            case OpCode::Sqrt:
                stack.push_back(std::sqrt(pop_one(stack, "sqrt")));
                break;
            case OpCode::Exp:
                stack.push_back(std::exp(pop_one(stack, "exp")));
                break;
            case OpCode::Log:
                stack.push_back(std::log(pop_one(stack, "log")));
                break;
            case OpCode::Sin:
                stack.push_back(std::sin(pop_one(stack, "sin")));
                break;
            case OpCode::Cos:
                stack.push_back(std::cos(pop_one(stack, "cos")));
                break;
            case OpCode::Abs:
                stack.push_back(std::abs(pop_one(stack, "abs")));
                break;
            case OpCode::Add:
            case OpCode::Subtract:
            case OpCode::Multiply:
            case OpCode::Divide:
            case OpCode::Power: {
                const double rhs = pop_one(stack, "binary operator");
                const double lhs = pop_one(stack, "binary operator");
                if (ins.op == OpCode::Add) {
                    stack.push_back(lhs + rhs);
                } else if (ins.op == OpCode::Subtract) {
                    stack.push_back(lhs - rhs);
                } else if (ins.op == OpCode::Multiply) {
                    stack.push_back(lhs * rhs);
                } else if (ins.op == OpCode::Divide) {
                    stack.push_back(lhs / rhs);
                } else {
                    stack.push_back(std::pow(lhs, rhs));
                }
                break;
            }
        }
        ensure_finite(stack.back(), source_);
    }

    if (stack.size() != 1) {
        throw std::invalid_argument("Expression evaluation ended with invalid stack state.");
    }
    ensure_finite(stack.back(), source_);
    return stack.back();
}

CompiledExpression compile_expression(const std::string& text) {
    return CompiledExpression(compile_to_rpn(tokenize(text)), text);
}

}  // namespace wop::config
