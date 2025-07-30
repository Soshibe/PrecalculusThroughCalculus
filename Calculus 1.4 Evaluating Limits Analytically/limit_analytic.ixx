export module limit_analytic;

import symbolic_function_parser;
import symbolic_function_differentiation;
import <memory>;
import <optional>;
import <stdexcept>;
import <cmath>;
import <string>;
import <limits>; // Required for std::numeric_limits

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
namespace limit_analytic {

    // A small epsilon for floating point comparisons, especially for checking zero.
    // It's important to use a consistent epsilon. This value should be carefully
    // chosen based on the expected precision needs.
    constexpr double EPSILON_LIMIT = 1e-9;

    // Evaluate f(x) by substituting x = val, returning optional<double> if evaluation succeeds.
    // Can now return +/- infinity, but nullopt for true indeterminate forms (0/0, inf/inf)
    // or unhandled symbolic parts.
    std::optional<double> try_substitute(NodePtr f, const std::string& var, double val);

    // Main analytic limit evaluator
    // Returns std::optional<NodePtr> representing the symbolic limit, or nullopt if not evaluable
    std::optional<NodePtr> analytic_limit(NodePtr f, const std::string& var, double val, int recursion_depth = 0);

    // Implementation details below:

    std::optional<double> try_substitute(NodePtr f, const std::string& var, double val) {
        if (!f) return std::nullopt;
        switch (f->type) {
        case FuncType::Constant:
            return f->value;
        case FuncType::Variable:
            if (f->name == var) return val;
            else return std::nullopt; // Variable other than var: cannot substitute to a numeric value
        case FuncType::Add: {
            auto left = try_substitute(f->children[0], var, val);
            auto right = try_substitute(f->children[1], var, val);
            if (left && right) return *left + *right;
            return std::nullopt;
        }
        case FuncType::Sub: {
            auto left = try_substitute(f->children[0], var, val);
            auto right = try_substitute(f->children[1], var, val);
            if (left && right) return *left - *right;
            return std::nullopt;
        }
        case FuncType::Mul: {
            auto left = try_substitute(f->children[0], var, val);
            auto right = try_substitute(f->children[1], var, val);
            if (left && right) return (*left) * (*right);
            return std::nullopt;
        }
        case FuncType::Div: {
            auto left_opt = try_substitute(f->children[0], var, val);
            auto right_opt = try_substitute(f->children[1], var, val);
            if (left_opt && right_opt) {
                double left = *left_opt;
                double right = *right_opt;

                // Check for division by zero
                if (std::fabs(right) < EPSILON_LIMIT) { // Denominator is effectively zero
                    if (std::fabs(left) < EPSILON_LIMIT) { // 0/0 form -> indeterminate
                        return std::nullopt; // Signal indeterminate form for L'Hopital
                    }
                    else { // X/0 form (X != 0) -> infinite limit
                        // Determine the sign of infinity. Divide by a tiny positive/negative value
                        // to get the correct sign, then apply to infinity.
                        double sign = (left > 0 && right > 0) || (left < 0 && right < 0) ? 1.0 : -1.0; // Simple sign based on numerator/denominator
                        // A more robust way to get sign, especially if 'right' is tiny and already signed zero
                        // For example: std::copysign(1.0, left) / std::copysign(0.0, right)
                        // Or simply std::copysign(std::numeric_limits<double>::infinity(), left / right);
                        return std::copysign(std::numeric_limits<double>::infinity(), left / right);
                    }
                }
                // Check for infinity/infinity form
                if (std::isinf(left) && std::isinf(right)) {
                    return std::nullopt; // Signal indeterminate form for L'Hopital
                }
                return left / right;
            }
            return std::nullopt; // Not enough numeric info to substitute
        }
        case FuncType::Pow: {
            auto base = try_substitute(f->children[0], var, val);
            auto exp = try_substitute(f->children[1], var, val);
            if (base && exp) {
                // Handle domain errors (e.g., negative base with fractional exponent)
                if (*base < 0.0 && std::floor(*exp) != *exp) return std::nullopt;
                // Handle 0^0, Inf^0, 1^Inf, Inf^0 etc. which are indeterminate forms
                if (std::fabs(*base) < EPSILON_LIMIT && std::fabs(*exp) < EPSILON_LIMIT) return std::nullopt; // 0^0
                if (std::isinf(*base) && std::fabs(*exp) < EPSILON_LIMIT) return std::nullopt; // Inf^0
                if (std::fabs(*base - 1.0) < EPSILON_LIMIT && std::isinf(*exp)) return std::nullopt; // 1^Inf

                return std::pow(*base, *exp);
            }
            return std::nullopt;
        }
        case FuncType::Neg: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) return -(*val_sub);
            return std::nullopt;
        }
        case FuncType::Sin: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) return std::sin(*val_sub);
            return std::nullopt;
        }
        case FuncType::Cos: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) return std::cos(*val_sub);
            return std::nullopt;
        }
        case FuncType::Tan: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) {
                // Check for values where tan is undefined (pi/2 + n*pi)
                double mod_val = std::fmod(std::fabs(*val_sub), M_PI);
                if (std::fabs(mod_val - M_PI_2) < EPSILON_LIMIT || std::fabs(mod_val - 3 * M_PI_2) < EPSILON_LIMIT) {
                    return std::nullopt; // Indicates undefined, leads to infinite limit. Let it fall through.
                }
                return std::tan(*val_sub);
            }
            return std::nullopt;
        }
        case FuncType::Ln: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub && *val_sub > EPSILON_LIMIT) return std::log(*val_sub);
            // If val_sub is 0 or negative, ln is undefined or -infinity.
            if (val_sub && *val_sub >= -EPSILON_LIMIT && *val_sub <= EPSILON_LIMIT) return -std::numeric_limits<double>::infinity(); // ln(0) = -inf
            return std::nullopt; // NaN or undefined
        }
        case FuncType::Log: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub && *val_sub > EPSILON_LIMIT) return std::log10(*val_sub);
            if (val_sub && *val_sub >= -EPSILON_LIMIT && *val_sub <= EPSILON_LIMIT) return -std::numeric_limits<double>::infinity(); // log10(0) = -inf
            return std::nullopt;
        }
        case FuncType::Abs: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) return std::abs(*val_sub);
            return std::nullopt;
        }
        case FuncType::Exp: {
            auto val_sub = try_substitute(f->children[0], var, val);
            if (val_sub) return std::exp(*val_sub);
            return std::nullopt;
        }
                          // For AddSubFunc, SubAddFunc, UnknownFunc etc., direct substitution usually isn't possible
                          // unless they simplify to a known form first.
        default:
            return std::nullopt;
        }
    }

    std::optional<NodePtr> analytic_limit(NodePtr f, const std::string& var, double val, int recursion_depth) {
        constexpr int MAX_RECURSION = 10; // To prevent infinite recursion with L'Hopital or simplification
        if (recursion_depth > MAX_RECURSION)
            throw std::runtime_error("Max recursion depth exceeded in analytic_limit.");

        if (!f) return std::nullopt;

        // 1. Try direct substitution first
        // try_substitute now returns +/- Inf for X/0, and nullopt for 0/0 or Inf/Inf or NaN cases.
        auto sub_val = try_substitute(f, var, val);
        if (sub_val.has_value()) {
            // If substitution yielded a finite number or +/- infinity, that's our limit.
            // std::isnan(*sub_val) is the specific case try_substitute returns nullopt for
            // 0/0, Inf/Inf and other truly indeterminate forms that need more analysis.
            // If try_substitute returns a NaN directly, it means it's an unresolvable
            // domain error (e.g. sqrt(-1)) rather than an indeterminate form.
            if (!std::isnan(*sub_val)) {
                return make_const(*sub_val);
            }
            // If it's NaN from try_substitute, we fall through, hoping L'Hopital or simplification helps.
            // However, it's more likely to remain unresolvable.
        }

        // 2. Apply L'Hôpital's Rule for 0/0 or ∞/∞ forms
        // This only applies if `f` is a division node.
        if (f->type == FuncType::Div) {
            auto numerator = f->children[0];
            auto denominator = f->children[1];

            // Evaluate numerator and denominator separately at the limit point.
            auto num_at_val_opt = try_substitute(numerator, var, val);
            auto den_at_val_opt = try_substitute(denominator, var, val);

            bool is_0_0 = false;
            bool is_inf_inf = false;

            if (num_at_val_opt && den_at_val_opt) {
                double num_val = *num_at_val_opt;
                double den_val = *den_at_val_opt;

                is_0_0 = (std::fabs(num_val) < EPSILON_LIMIT && std::fabs(den_val) < EPSILON_LIMIT);
                is_inf_inf = (std::isinf(num_val) && std::isinf(den_val));
            }

            if (is_0_0 || is_inf_inf) {
                auto d_num = derivative(numerator, var);
                auto d_den = derivative(denominator, var);

                // L'Hopital's applies if derivatives can be computed.
                if (d_num && d_den) {
                    // Form new function d_num/d_den and try limit recursively
                    auto new_f = make_op("/", d_num, d_den);
                    return analytic_limit(new_f, var, val, recursion_depth + 1);
                }
                // If derivatives can't be computed symbolically, or are null.
                return std::nullopt;
            }
        }

        // 3. Try to simplify the expression and retry.
        // This is crucial for expressions that might simplify to a solvable form
        // before or after direct substitution or L'Hopital (e.g., (x^2-1)/(x-1) -> x+1).
        NodePtr simplified = simplify(f);
        if (!are_equal(simplified, f)) { // Check if simplification actually changed the node
            // Recurse with the simplified form, incrementing depth to avoid infinite loops.
            return analytic_limit(simplified, var, val, recursion_depth + 1);
        }

        // 4. If all fails, return nullopt.
        // This covers limits that genuinely don't exist (e.g., oscillating functions),
        // expressions with unhandled function types or other variables,
        // and indeterminate forms not covered by L'Hopital (e.g., `inf - inf`, `0 * inf`).
        return std::nullopt;
    }

} // namespace limit_analytic

export using limit_analytic::analytic_limit;