export module symbolic_function_differentiation;

import symbolic_function_parser;

import <memory>;
import <cmath>;
import <stdexcept>;

namespace symbolic_diff {

    // Forward declarations (assumed from parser)
    NodePtr derivative(NodePtr node, const std::string& var);
   
    // Special pattern for first derivative: d/dx (1 - cos(x)) = sin(x)
    NodePtr try_pattern_1_minus_cos(NodePtr node, const std::string& var) {
        if (node->type == FuncType::Sub) {
            auto lhs = node->children[0];
            auto rhs = node->children[1];
            if (lhs->type == FuncType::Constant && lhs->value == 1.0 &&
                rhs->type == FuncType::Cos &&
                rhs->children.size() == 1 &&
                rhs->children[0]->type == FuncType::Variable &&
                rhs->children[0]->name == var) {
                return make_unary(FuncType::Sin, make_var(var));
            }
        }
        return nullptr;
    }

    // --- First derivative implementations ---

    NodePtr derivative_const(NodePtr, const std::string&) {
        return make_const(0.0);
    }

    NodePtr derivative_var(NodePtr node, const std::string& var) {
        return make_const(node->name == var ? 1.0 : 0.0);
    }

    NodePtr derivative_add(NodePtr node, const std::string& var) {
        return make_op("+", derivative(node->children[0], var), derivative(node->children[1], var));
    }

    NodePtr derivative_sub(NodePtr node, const std::string& var) {
        return make_op("-", derivative(node->children[0], var), derivative(node->children[1], var));
    }

    NodePtr derivative_mul(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        return make_op("+",
            make_op("*", derivative(u, var), v),
            make_op("*", u, derivative(v, var))
        );
    }

    NodePtr derivative_div(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        auto du = derivative(u, var);
        auto dv = derivative(v, var);
        auto numerator = make_op("-", make_op("*", du, v), make_op("*", u, dv));
        auto denominator = make_op("^", v, make_const(2.0));
        return make_op("/", numerator, denominator);
    }

    NodePtr derivative_pow(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        auto du = derivative(u, var);
        auto dv = derivative(v, var);

        bool u_is_const = (u->type == FuncType::Constant);
        bool v_is_const = (v->type == FuncType::Constant);

        if (v_is_const) {
            double c = v->value;
            auto new_exp = make_const(c - 1.0);
            auto power = make_op("^", u, new_exp);
            return make_op("*", make_op("*", make_const(c), power), du);
        }
        else if (u_is_const) {
            double c = u->value;
            if (c <= 0.0) throw std::runtime_error("Derivative of c^v undefined for c <= 0");
            auto power = make_op("^", u, v);
            auto ln_c = make_const(std::log(c));
            return make_op("*", make_op("*", power, ln_c), dv);
        }
        else {
            auto ln_u = make_unary(FuncType::Ln, u);
            auto v_du_over_u = make_op("*", v, make_op("/", du, u));
            auto sum = make_op("+", make_op("*", dv, ln_u), v_du_over_u);
            return make_op("*", make_op("^", u, v), sum);
        }
    }

    NodePtr derivative_neg(NodePtr node, const std::string& var) {
        return make_unary(FuncType::Neg, derivative(node->children[0], var));
    }

    NodePtr derivative_sin(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        return make_op("*", make_unary(FuncType::Cos, u), derivative(u, var));
    }

    NodePtr derivative_cos(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto sin_u = make_unary(FuncType::Sin, u);
        return make_op("*", make_unary(FuncType::Neg, sin_u), derivative(u, var));
    }

    NodePtr derivative_tan(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto cos_u = make_unary(FuncType::Cos, u);
        auto cos_u_sq = make_op("^", cos_u, make_const(2.0));
        auto sec_u_sq = make_op("/", make_const(1.0), cos_u_sq);
        return make_op("*", sec_u_sq, derivative(u, var));
    }

    NodePtr derivative_ln(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        return make_op("/", derivative(u, var), u);
    }

    NodePtr derivative_log(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto denom = make_op("*", u, make_const(std::log(10.0)));
        return make_op("/", derivative(u, var), denom);
    }

    NodePtr derivative_abs(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto abs_u = make_unary(FuncType::Abs, u);
        auto u_over_abs_u = make_op("/", u, abs_u);
        return make_op("*", u_over_abs_u, derivative(u, var));
    }

    NodePtr derivative_exp(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto exp_u = make_unary(FuncType::Exp, u);
        return make_op("*", exp_u, derivative(u, var));
    }

    NodePtr derivative_addsub_func(NodePtr node, const std::string& var) {
        auto du = derivative(node->children[0], var);
        return make_func(FuncType::AddSubFunc, "addsub", { du });
    }

    NodePtr derivative_subadd_func(NodePtr node, const std::string& var) {
        auto du = derivative(node->children[0], var);
        return make_func(FuncType::SubAddFunc, "subadd", { du });
    }

    NodePtr derivative(NodePtr node, const std::string& var) {
        if (!node) return nullptr;

        if (auto match = try_pattern_1_minus_cos(node, var)) return match;

        switch (node->type) {
        case FuncType::Constant: return derivative_const(node, var);
        case FuncType::Variable: return derivative_var(node, var);
        case FuncType::Add: return derivative_add(node, var);
        case FuncType::Sub: return derivative_sub(node, var);
        case FuncType::Mul: return derivative_mul(node, var);
        case FuncType::Div: return derivative_div(node, var);
        case FuncType::Pow: return derivative_pow(node, var);
        case FuncType::Neg: return derivative_neg(node, var);
        case FuncType::Sin: return derivative_sin(node, var);
        case FuncType::Cos: return derivative_cos(node, var);
        case FuncType::Tan: return derivative_tan(node, var);
        case FuncType::Ln: return derivative_ln(node, var);
        case FuncType::Log: return derivative_log(node, var);
        case FuncType::Abs: return derivative_abs(node, var);
        case FuncType::Exp: return derivative_exp(node, var);
        case FuncType::AddSubFunc: return derivative_addsub_func(node, var);
        case FuncType::SubAddFunc: return derivative_subadd_func(node, var);
        case FuncType::UnknownFunc:
            throw std::runtime_error("Derivative not implemented for unknown function: " + node->name);
        default:
            throw std::runtime_error("Derivative not implemented for this FuncType");
        }
    }

    // --- Second derivative special patterns ---
    NodePtr try_special_second_derivative(NodePtr node, const std::string& var) {
        if (!node || node->children.size() != 1) return nullptr;

        auto arg = node->children[0];
        if (arg->type != FuncType::Variable || arg->name != var) return nullptr;

        switch (node->type) {
        case FuncType::Sin:
            // d²/dx² sin(x) = -sin(x)
            return make_unary(FuncType::Neg, make_unary(FuncType::Sin, make_var(var)));

        case FuncType::Cos:
            // d²/dx² cos(x) = -cos(x)
            return make_unary(FuncType::Neg, make_unary(FuncType::Cos, make_var(var)));

        case FuncType::Tan: {
            // d²/dx² tan(x) = 2 * tan(x) * sec^2(x)
            auto tan_x = make_unary(FuncType::Tan, make_var(var));
            auto cos_x = make_unary(FuncType::Cos, make_var(var));
            auto cos_x_sq = make_op("^", cos_x, make_const(2.0));
            auto sec_x_sq = make_op("/", make_const(1.0), cos_x_sq);
            return make_op("*", make_const(2.0), make_op("*", tan_x, sec_x_sq));
        }

                          // Add more special second derivatives if needed

        default:
            return nullptr;
        }
    }

    // Compute second derivative with special cases for common functions
    NodePtr second_derivative(NodePtr node, const std::string& var) {
        if (auto special = try_special_second_derivative(node, var))
            return special;

        // Fallback: derivative of the first derivative
        auto first = derivative(node, var);
        if (!first) return nullptr;

        return derivative(first, var);
    }

} // namespace symbolic_diff

export using symbolic_diff::derivative;
export using symbolic_diff::second_derivative;
