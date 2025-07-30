export module symbolic_function_differentiation;

import symbolic_function_parser;

import <memory>;
import <cmath>;
import <stdexcept>;

namespace symbolic_diff {

    // Forward declarations (already mostly defined by 'using namespace')
    NodePtr derivative(NodePtr node, const std::string& var);
    // NodePtr make_const(double val); // Defined in parser module, brought by 'using namespace'
    // NodePtr make_var(const std::string& name); // Defined in parser module, brought by 'using namespace'
    // NodePtr make_op(const std::string& op, NodePtr lhs, NodePtr rhs); // Defined in parser module, brought by 'using namespace'
    // NodePtr make_unary(FuncType ft, NodePtr arg); // Defined in parser module, brought by 'using namespace'
    // NodePtr simplify(NodePtr node); // Defined in parser module, brought by 'using namespace'

    // d/dx (constant) = 0
    NodePtr derivative_const(NodePtr node, const std::string&) {
        return make_const(0.0);
    }

    // d/dx (variable) = 1 if variable matches, else 0
    NodePtr derivative_var(NodePtr node, const std::string& var) {
        return make_const(node->name == var ? 1.0 : 0.0);
    }

    // Sum rule: d/dx (u + v) = du/dx + dv/dx
    NodePtr derivative_add(NodePtr node, const std::string& var) {
        return make_op("+", derivative(node->children[0], var), derivative(node->children[1], var));
    }

    // Difference rule: d/dx (u - v) = du/dx - dv/dx
    NodePtr derivative_sub(NodePtr node, const std::string& var) {
        return make_op("-", derivative(node->children[0], var), derivative(node->children[1], var));
    }

    // Product rule: d/dx (u * v) = u' * v + u * v'
    NodePtr derivative_mul(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        auto du = derivative(u, var);
        auto dv = derivative(v, var);
        return make_op("+", make_op("*", du, v), make_op("*", u, dv));
    }

    // Quotient rule: d/dx (u / v) = (u' * v - u * v') / v^2
    NodePtr derivative_div(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        auto du = derivative(u, var);
        auto dv = derivative(v, var);
        auto numerator = make_op("-", make_op("*", du, v), make_op("*", u, dv));
        auto denominator = make_op("^", v, make_const(2.0));
        return make_op("/", numerator, denominator);
    }

    // Power rule and general power differentiation:
    // For f(x) = u(x)^v(x), derivative is:
    // f' = u^v * [v' * ln(u) + v * u'/u]
    NodePtr derivative_pow(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto v = node->children[1];
        auto du = derivative(u, var);
        auto dv = derivative(v, var);

        bool u_is_const = (u->type == FuncType::Constant);
        bool v_is_const = (v->type == FuncType::Constant);

        if (v_is_const) {
            // Simplified power rule: d/dx u^c = c * u^(c-1) * u'
            double c = v->value;
            auto new_exp = make_const(c - 1.0);
            auto power = make_op("^", u, new_exp);
            return make_op("*", make_op("*", make_const(c), power), du);
        }
        else if (u_is_const) {
            // d/dx c^v = c^v * ln(c) * v'
            double c = u->value;
            if (c <= 0.0) throw std::runtime_error("Derivative of c^v undefined for c <= 0");
            auto power = make_op("^", u, v);
            auto ln_c = make_const(std::log(c));
            return make_op("*", make_op("*", power, ln_c), dv);
        }
        else {
            // General case: u^v * (v' * ln(u) + v * u'/u)
            auto ln_u = make_unary(FuncType::Ln, u);
            auto v_du_over_u = make_op("*", v, make_op("/", du, u));
            auto sum = make_op("+", make_op("*", dv, ln_u), v_du_over_u);
            return make_op("*", make_op("^", u, v), sum);
        }
    }

    // Derivative of unary minus: d/dx (-u) = -u'
    NodePtr derivative_neg(NodePtr node, const std::string& var) {
        return make_unary(FuncType::Neg, derivative(node->children[0], var));
    }

    // Derivatives of elementary functions:
    // d/dx sin(u) = cos(u) * u'
    NodePtr derivative_sin(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        return make_op("*", make_unary(FuncType::Cos, u), du);
    }

    // d/dx cos(u) = -sin(u) * u'
    NodePtr derivative_cos(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        auto sin_u = make_unary(FuncType::Sin, u);
        return make_op("*", make_unary(FuncType::Neg, sin_u), du);
    }

    // d/dx tan(u) = sec^2(u) * u'
    // sec^2(u) = 1 / cos^2(u)
    NodePtr derivative_tan(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        auto cos_u = make_unary(FuncType::Cos, u);
        auto cos_u_sq = make_op("^", cos_u, make_const(2.0));
        auto sec_u_sq = make_op("/", make_const(1.0), cos_u_sq);
        return make_op("*", sec_u_sq, du);
    }

    // d/dx ln(u) = u'/u
    NodePtr derivative_ln(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        return make_op("/", du, u);
    }

    // d/dx log(u) = u' / (u * ln(10))
    NodePtr derivative_log(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        auto denom = make_op("*", u, make_const(std::log(10.0)));
        return make_op("/", du, denom);
    }

    // d/dx abs(u) = (u / abs(u)) * u' , undefined at u=0 (leave as is)
    NodePtr derivative_abs(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        auto abs_u = make_unary(FuncType::Abs, u);
        auto u_over_abs_u = make_op("/", u, abs_u);
        return make_op("*", u_over_abs_u, du);
    }

    // d/dx exp(u) = exp(u) * u'
    NodePtr derivative_exp(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        auto exp_u = make_unary(FuncType::Exp, u);
        return make_op("*", exp_u, du);
    }

    // New: d/dx addsub(u) = addsub(u')
    NodePtr derivative_addsub_func(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        return make_func(FuncType::AddSubFunc, "addsub", { du });
    }

    // New: d/dx subadd(u) = subadd(u')
    NodePtr derivative_subadd_func(NodePtr node, const std::string& var) {
        auto u = node->children[0];
        auto du = derivative(u, var);
        return make_func(FuncType::SubAddFunc, "subadd", { du });
    }


    // Main derivative dispatcher
    NodePtr derivative(NodePtr node, const std::string& var) {
        if (!node) return nullptr;

        switch (node->type) {
        case FuncType::Constant:
            return derivative_const(node, var);
        case FuncType::Variable:
            return derivative_var(node, var);
        case FuncType::Add:
            return derivative_add(node, var);
        case FuncType::Sub:
            return derivative_sub(node, var);
        case FuncType::Mul:
            return derivative_mul(node, var);
        case FuncType::Div:
            return derivative_div(node, var);
        case FuncType::Pow:
            return derivative_pow(node, var);
        case FuncType::Neg:
            return derivative_neg(node, var);
        case FuncType::Sin:
            return derivative_sin(node, var);
        case FuncType::Cos:
            return derivative_cos(node, var);
        case FuncType::Tan:
            return derivative_tan(node, var);
        case FuncType::Ln:
            return derivative_ln(node, var);
        case FuncType::Log:
            return derivative_log(node, var);
        case FuncType::Abs:
            return derivative_abs(node, var);
        case FuncType::Exp:
            return derivative_exp(node, var);
        case FuncType::AddSubFunc: // Handle new function type
            return derivative_addsub_func(node, var);
        case FuncType::SubAddFunc: // Handle new function type
            return derivative_subadd_func(node, var);
        case FuncType::UnknownFunc:
            // For unknown functions, we assume they are not dependent on 'var'
            // unless we have specific rules. For a general symbolic diff,
            // treating them as constants with respect to 'var' is often
            // the default, or throwing an error if strict.
            // For now, let's treat them as constants (derivative is 0)
            // or if they represent a function f(var), we'd need f'(var).
            // A safer default is to return an error or a symbolic derivative f'(x)
            // For simplicity, for an unknown function name, we assume it's
            // not directly dependent on `var` if it's not explicitly passed
            // as an argument. If it's `f(x)`, then d/dx f(x) is f'(x)
            // If it's `f(y)`, then d/dx f(y) is 0.
            // Given it's a unary function of one argument, the chain rule
            // applies to the argument: d/dx(f(u)) = f'(u) * u'.
            // However, we don't know f'. So, for now, we'll throw an error.
            throw std::runtime_error("Derivative not implemented for unknown function type or custom function name: " + node->name);
        default:
            throw std::runtime_error("Derivative not implemented for this FuncType");
        }
    }

} // namespace symbolic_diff

export using symbolic_diff::derivative;