export module limit_analytic;

import symbolic_function_parser;
import symbolic_function_differentiation;
import symbolic_polynomial_analysis;
import <memory>;
import <optional>;
import <stdexcept>; // For std::runtime_error
import <cmath>;     // For std::fabs, std::isinf, std::copysign, std::sin, std::cos, std::tan, std::log, std::log10, std::exp, std::pow, std::fmod
import <string>;    // For std::string, std::to_string
import <limits>;    // Required for std::numeric_limits
import <variant>;   // For std::variant
import <iostream>;  // For std::cerr

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

namespace limit_analytic {
    export std::string funcTypeToString(FuncType type) {
        switch (type) {
        case FuncType::Constant: return "Constant";
        case FuncType::Variable: return "Variable";
        case FuncType::Add: return "+";
        case FuncType::Sub: return "-";
        case FuncType::Mul: return "*";
        case FuncType::Div: return "/";
        case FuncType::Pow: return "^";
        case FuncType::Neg: return "-";
        case FuncType::Sin: return "sin";
        case FuncType::Cos: return "cos";
        case FuncType::Tan: return "tan";
        case FuncType::Log: return "log";
        case FuncType::Ln: return "ln";
        case FuncType::Abs: return "abs";
        case FuncType::Exp: return "exp";
        case FuncType::AddSubFunc: return "AddSubFunc";
        case FuncType::SubAddFunc: return "SubAddFunc";
        case FuncType::UnknownFunc: return "Unknown Function Type";
        default: return "Invalid FuncType Value";
        }
    }

    constexpr double EPSILON_LIMIT = 1e-9;
    constexpr double infinity = std::numeric_limits<double>::infinity();

    export enum class LimitPointType {
        Finite,
        PositiveInfinity,
        NegativeInfinity
    };

    export enum class LimitDirection {
        Both,
        FromRight,
        FromLeft
    };

    export struct LimitPoint {
        LimitPointType type;
        double value;
        std::optional<LimitDirection> direction;

        static LimitPoint finite(double v, std::optional<LimitDirection> dir = std::nullopt) {
            if (std::fabs(v) < EPSILON_LIMIT && !dir.has_value()) {
                dir = LimitDirection::Both;
            }
            return { LimitPointType::Finite, v, dir };
        }
        static LimitPoint positive_infinity() { return { LimitPointType::PositiveInfinity, infinity, std::nullopt }; }
        static LimitPoint negative_infinity() { return { LimitPointType::NegativeInfinity, -infinity, std::nullopt }; }
    };

    export LimitPoint parse_limit_point_string(const std::string& val_str) {
        if (val_str == "inf" || val_str == "infinity" || val_str == "+inf" || val_str == "+infinity") {
            return LimitPoint::positive_infinity();
        }
        if (val_str == "-inf" || val_str == "-infinity") {
            return LimitPoint::negative_infinity();
        }
        if (val_str == "0+") {
            return LimitPoint::finite(0.0, LimitDirection::FromRight);
        }
        if (val_str == "0-") {
            return LimitPoint::finite(0.0, LimitDirection::FromLeft);
        }
        try {
            size_t pos;
            double val = std::stod(val_str, &pos);
            if (pos == val_str.length()) {
                return LimitPoint::finite(val);
            }
        }
        catch (const std::invalid_argument&) {
        }
        catch (const std::out_of_range&) {
        }
        throw std::runtime_error("Invalid limit point string: " + val_str + ". Expected a number, 'inf', '-inf', '0+', or '0-'.");
    }

    std::optional<double> try_substitute(NodePtr f, const std::string& var, LimitPoint val_point);

    export std::optional<NodePtr> analytic_limit(NodePtr f, const std::string& var, LimitPoint val_point, int recursion_depth = 0);

    std::optional<double> evaluate_with_perturbation(NodePtr node, const std::string& var, LimitPoint val_point) {
        if (val_point.type != LimitPointType::Finite || !val_point.direction.has_value() || val_point.direction.value() == LimitDirection::Both) {
            return try_substitute(node, var, val_point);
        }

        double perturbation_offset;
        if (val_point.direction.value() == LimitDirection::FromRight) {
            perturbation_offset = EPSILON_LIMIT;
        }
        else {
            perturbation_offset = -EPSILON_LIMIT;
        }

        LimitPoint perturbed_limit_point = LimitPoint::finite(val_point.value + perturbation_offset, val_point.direction);

        return try_substitute(node, var, perturbed_limit_point);
    }

    int count_nodes(NodePtr node) {
        if (!node) return 0;
        int count = 1;
        for (const auto& child : node->children) {
            count += count_nodes(child);
        }
        return count;
    }

    std::string to_string(NodePtr node) {
        if (!node) return "null";
        switch (node->type) {
        case FuncType::Constant: return std::to_string(node->value);
        case FuncType::Variable: return node->name;
        case FuncType::Neg: return "-" + to_string(node->children[0]);
        case FuncType::Add:
        case FuncType::Sub:
        case FuncType::Mul:
        case FuncType::Div:
        case FuncType::Pow:
            return "(" + to_string(node->children[0]) + " " + limit_analytic::funcTypeToString(node->type) + " " + to_string(node->children[1]) + ")";
        case FuncType::Sin:
        case FuncType::Cos:
        case FuncType::Tan:
        case FuncType::Log:
        case FuncType::Ln:
        case FuncType::Abs:
        case FuncType::Exp:
        case FuncType::AddSubFunc:
        case FuncType::SubAddFunc:
            if (!node->children.empty()) {
                return limit_analytic::funcTypeToString(node->type) + "(" + to_string(node->children[0]) + ")";
            }
            return limit_analytic::funcTypeToString(node->type) + "()";
        default: return "UNKNOWN_NODE_TYPE";
        }
    }


    std::optional<double> try_substitute(NodePtr f, const std::string& var, LimitPoint val_point) {
        if (!f) return std::nullopt;

        double finite_val = val_point.value;
        bool is_inf_limit = (val_point.type == LimitPointType::PositiveInfinity || val_point.type == LimitPointType::NegativeInfinity);
        bool is_pos_inf_limit = (val_point.type == LimitPointType::PositiveInfinity);

        switch (f->type) {
        case FuncType::Constant:
            return f->value;
        case FuncType::Variable:
            if (f->name == var) {
                if (is_inf_limit) return val_point.value;
                else return finite_val;
            }
            else return std::nullopt;
        case FuncType::Add: {
            auto left_opt = try_substitute(f->children[0], var, val_point);
            auto right_opt = try_substitute(f->children[1], var, val_point);

            if (!left_opt || !right_opt) return std::nullopt;

            double left = *left_opt;
            double right = *right_opt;

            if (std::isinf(left) && std::isinf(right)) {
                if (std::copysign(1.0, left) != std::copysign(1.0, right)) {
                    return std::nullopt;
                }
                return left;
            }
            if (std::isinf(left)) return left;
            if (std::isinf(right)) return right;

            return left + right;
        }
        case FuncType::Sub: {
            auto left_opt = try_substitute(f->children[0], var, val_point);
            auto right_opt = try_substitute(f->children[1], var, val_point);

            if (!left_opt || !right_opt) return std::nullopt;

            double left = *left_opt;
            double right = *right_opt;

            if (std::isinf(left) && std::isinf(right)) {
                if (std::copysign(1.0, left) == std::copysign(1.0, right)) {
                    return std::nullopt;
                }
                return left;
            }
            if (std::isinf(left)) return left;
            if (std::isinf(right)) return -right;

            return left - right;
        }
        case FuncType::Mul: {
            auto left_opt = try_substitute(f->children[0], var, val_point);
            auto right_opt = try_substitute(f->children[1], var, val_point);

            if (!left_opt || !right_opt) return std::nullopt;

            double left = *left_opt;
            double right = *right_opt;

            if (std::isinf(left) && std::isinf(right)) {
                return std::copysign(infinity, left * right);
            }
            if (std::isinf(left) && std::fabs(right) < EPSILON_LIMIT) return std::nullopt;
            if (std::isinf(right) && std::fabs(left) < EPSILON_LIMIT) return std::nullopt;
            if (std::isinf(left)) return std::copysign(infinity, left * right);
            if (std::isinf(right)) return std::copysign(infinity, left * right);

            return left * right;
        }
                          // ... (rest of the try_substitute function remains the same until FuncType::Div)

        case FuncType::Div: {
            auto left_opt = try_substitute(f->children[0], var, val_point);
            auto right_opt = try_substitute(f->children[1], var, val_point);

            if (!left_opt || !right_opt) return std::nullopt;

            double left = *left_opt;
            double right = *right_opt;

            bool is_0_0 = (std::fabs(left) < EPSILON_LIMIT && std::fabs(right) < EPSILON_LIMIT);
            bool is_inf_inf = (std::isinf(left) && std::isinf(right));

            if (is_0_0 || is_inf_inf) {
                // This is an indeterminate form, cannot resolve by direct substitution.
                return std::nullopt;
            }

            if (std::fabs(right) < EPSILON_LIMIT) { // Denominator is zero (or very close to zero)
                // Numerator is non-zero, denominator is zero: results in infinity.
                // We need to determine the sign of infinity using a perturbed value for the denominator.
                auto perturbed_den_opt = evaluate_with_perturbation(f->children[1], var, val_point);

                if (perturbed_den_opt.has_value() && std::fabs(*perturbed_den_opt) > EPSILON_LIMIT) {
                    // If perturbation helps to make the denominator non-zero, calculate the signed infinity.
                    // This handles cases like 1/x as x->0+ or x->0-.
                    return std::copysign(infinity, left * (*perturbed_den_opt));
                }
                else if (perturbed_den_opt.has_value() && std::fabs(*perturbed_den_opt) < EPSILON_LIMIT) {
                    // If the denominator is still zero after perturbation, it's an unresolvable division by zero.
                    // This scenario should be rare if perturbation is effective.
                    std::cerr << "WARNING: Denominator still zero after perturbation in Div for finite/zero case. Returning nullopt.\n";
                    return std::nullopt; // Or throw an error/handle differently if this case should never happen.
                }
                else {
                    // Perturbation did not yield a value (e.g., if the perturbed expression is undefined).
                    return std::nullopt;
                }
            }

            // Standard division if none of the above special cases.
            if (!std::isinf(left) && !std::isinf(right)) {
                return left / right;
            }

            // Cases involving infinity where not indeterminate (e.g., finite/inf -> 0, inf/finite -> inf)
            if (!std::isinf(left) && std::isinf(right)) {
                return 0.0;
            }

            if (std::isinf(left) && !std::isinf(right)) {
                return std::copysign(infinity, left / right);
            }

            // Fallback for any other unhandled combinations involving infinity.
            return std::nullopt;
        }
                          // ... (rest of the try_substitute function)
        case FuncType::Pow: {
            auto base_opt = try_substitute(f->children[0], var, val_point);
            auto exp_opt = try_substitute(f->children[1], var, val_point);

            if (!base_opt || !exp_opt) return std::nullopt;

            double base = *base_opt;
            double exp = *exp_opt;

            if (std::fabs(base) < EPSILON_LIMIT && std::fabs(exp) < EPSILON_LIMIT) return std::nullopt;
            if (std::fabs(base - 1.0) < EPSILON_LIMIT && std::isinf(exp)) return std::nullopt;
            if (std::isinf(base) && std::fabs(exp) < EPSILON_LIMIT) return std::nullopt;

            if (std::isinf(exp)) {
                if (exp > 0) {
                    if (base > 1.0 + EPSILON_LIMIT) return infinity;
                    if (base > 0.0 && base < 1.0 - EPSILON_LIMIT) return 0.0;
                    if (base < -1.0 - EPSILON_LIMIT) return infinity;
                    if (base >= -1.0 && base <= -EPSILON_LIMIT) return 0.0;
                    if (std::fabs(base - 1.0) < EPSILON_LIMIT) return 1.0;
                    if (std::fabs(base + 1.0) < EPSILON_LIMIT) return std::nullopt;
                }
                else {
                    if (base > 1.0 + EPSILON_LIMIT) return 0.0;
                    if (base > 0.0 && base < 1.0 - EPSILON_LIMIT) return infinity;
                    if (base < -1.0 - EPSILON_LIMIT) return 0.0;
                    if (base >= -1.0 && base <= -EPSILON_LIMIT) return infinity;
                    if (std::fabs(base - 1.0) < EPSILON_LIMIT) return 1.0;
                    if (std::fabs(base + 1.0) < EPSILON_LIMIT) return std::nullopt;
                }
            }

            if (base < 0.0 && std::floor(exp) != exp) return std::nullopt;
            if (base == 0.0 && exp < 0.0) return infinity;

            return std::pow(base, exp);
        }
        case FuncType::Neg: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) return -(*val_sub);
            return std::nullopt;
        }
        case FuncType::Sin: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) return std::nullopt;
                return std::sin(*val_sub);
            }
            return std::nullopt;
        }
        case FuncType::Cos: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) return std::nullopt;
                return std::cos(*val_sub);
            }
            return std::nullopt;
        }
        case FuncType::Tan: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) return std::nullopt;
                double mod_val = std::fmod(std::fabs(*val_sub), M_PI);
                bool is_odd_pi_half = (std::fabs(mod_val - M_PI_2) < EPSILON_LIMIT);
                if (is_odd_pi_half) {
                    return std::nullopt;
                }
                return std::tan(*val_sub);
            }
            return std::nullopt;
        }
        case FuncType::Ln: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) {
                    if (*val_sub > 0) return infinity;
                    return std::nullopt;
                }
                if (std::fabs(*val_sub) < EPSILON_LIMIT) {
                    auto perturbed_arg_opt = evaluate_with_perturbation(f->children[0], var, val_point);

                    if (perturbed_arg_opt.has_value()) {
                        if (*perturbed_arg_opt > EPSILON_LIMIT) {
                            return -infinity;
                        }
                        else if (*perturbed_arg_opt < -EPSILON_LIMIT) {
                            return std::nullopt;
                        }
                    }
                    return std::nullopt;
                }
                if (*val_sub > EPSILON_LIMIT) return std::log(*val_sub);
                return std::nullopt;
            }
            return std::nullopt;
        }
        case FuncType::Log: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) {
                    if (*val_sub > 0) return infinity;
                    return std::nullopt;
                }
                if (std::fabs(*val_sub) < EPSILON_LIMIT) {
                    auto perturbed_arg_opt = evaluate_with_perturbation(f->children[0], var, val_point);

                    if (perturbed_arg_opt.has_value()) {
                        if (*perturbed_arg_opt > EPSILON_LIMIT) {
                            return -infinity;
                        }
                        else if (*perturbed_arg_opt < -EPSILON_LIMIT) {
                            return std::nullopt;
                        }
                    }
                    return std::nullopt;
                }
                if (*val_sub > EPSILON_LIMIT) return std::log10(*val_sub);
                return std::nullopt;
            }
            return std::nullopt;
        }
        case FuncType::Abs: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) return std::abs(*val_sub);
            return std::nullopt;
        }
        case FuncType::Exp: {
            auto val_sub = try_substitute(f->children[0], var, val_point);
            if (val_sub) {
                if (std::isinf(*val_sub)) {
                    if (*val_sub > 0) return infinity;
                    return 0.0;
                }
                return std::exp(*val_sub);
            }
            return std::nullopt;
        }
        case FuncType::AddSubFunc:
        case FuncType::SubAddFunc:
        case FuncType::UnknownFunc:
        default:
            return std::nullopt;
        }
    }

    std::optional<NodePtr> analytic_limit(NodePtr f, const std::string& var, LimitPoint val_point, int recursion_depth) {
        constexpr int MAX_RECURSION = 10;
        constexpr int MAX_NODE_COMPLEXITY = 500;

        if (recursion_depth > MAX_RECURSION) {
            std::cerr << "DEBUG: Max recursion depth (" << MAX_RECURSION << ") exceeded. Returning nullopt.\n";
            return std::nullopt;
        }

        if (!f) return std::nullopt;

        std::string expression_str = to_string(f);

        std::cerr << "DEBUG: Entering analytic_limit for expression: " << expression_str
            << ", as " << var << " -> " << (val_point.type == LimitPointType::Finite ? std::to_string(val_point.value) : (val_point.type == LimitPointType::PositiveInfinity ? "inf" : "-inf"))
            << (val_point.direction.has_value() && val_point.direction != LimitDirection::Both ? (val_point.direction == LimitDirection::FromRight ? "+" : "-") : "")
            << ", recursion_depth: " << recursion_depth << "\n";

        // 1. Always try to simplify the expression first.
        NodePtr simplified = simplify_recursive(f);
        if (!are_equal(simplified, f)) {
            std::cerr << "DEBUG: Simplified expression from " << to_string(f) << " to " << to_string(simplified) << ". Retrying limit evaluation with recursion_depth: " << recursion_depth + 1 << ".\n";
            // Recurse with the simplified expression, incrementing depth to prevent infinite loops if simplification loops
            return analytic_limit(simplified, var, val_point, recursion_depth + 1);
        }
        // If no simplification occurred, proceed with the original function 'f'.

        // 2. Try direct substitution
        auto sub_val = try_substitute(f, var, val_point);
        if (sub_val.has_value()) {
            if (!std::isnan(*sub_val)) {
                std::cerr << "DEBUG: Direct substitution successful. Value: " << *sub_val << "\n";
                return make_const(*sub_val);
            }
            std::cerr << "DEBUG: Direct substitution resulted in NaN (likely an indeterminate form that try_substitute didn't catch as nullopt explicitly). Attempting other methods.\n";
        }
        else {
            std::cerr << "DEBUG: Direct substitution returned nullopt (indeterminate form or unhandled symbolic element).\n";
        }

        // 3. Apply L'Hôpital's Rule for 0/0 or ∞/∞ forms
        if (f->type == FuncType::Div) {
            NodePtr numerator = f->children[0];
            NodePtr denominator = f->children[1];

            auto num_at_val_opt = try_substitute(numerator, var, val_point);
            auto den_at_val_opt = try_substitute(denominator, var, val_point);

            bool is_0_0 = false;
            bool is_inf_inf = false;

            if (num_at_val_opt && den_at_val_opt) {
                double num_val = *num_at_val_opt;
                double den_val = *den_at_val_opt;

                is_0_0 = (std::fabs(num_val) < EPSILON_LIMIT && std::fabs(den_val) < EPSILON_LIMIT);
                is_inf_inf = (std::isinf(num_val) && std::isinf(den_val));

                std::cerr << "DEBUG: L'Hopital check: num_val=" << num_val << ", den_val=" << den_val
                    << " (is_0_0=" << is_0_0 << ", is_inf_inf=" << is_inf_inf << ")\n";

                if (is_0_0 || is_inf_inf) {
                    NodePtr d_num = derivative(numerator, var);
                    NodePtr d_den = derivative(denominator, var);

                    if (!d_num || !d_den) {
                        std::cerr << "DEBUG: Could not differentiate for L'Hopital's. Returning nullopt.\n";
                        return std::nullopt;
                    }

                    if (recursion_depth == 0) {
                        std::cerr << "L'Hopital's Rule: Applying first derivative (1st application).\n";
                    }
                    else if (recursion_depth == 1) {
                        std::cerr << "L'Hopital's Rule: Applying second derivative (2nd application).\n";
                        int num_nodes_d_num = count_nodes(d_num);
                        int num_nodes_d_den = count_nodes(d_den);

                        if (num_nodes_d_num > MAX_NODE_COMPLEXITY || num_nodes_d_den > MAX_NODE_COMPLEXITY) {
                            std::cerr << "WARNING: Second L'Hopital application might lead to expression explosion.\n";
                            std::cerr << "        Nodes in d_num: " << num_nodes_d_num
                                << ", Nodes in d_den: " << num_nodes_d_den << " (Threshold: " << MAX_NODE_COMPLEXITY << ")\n";
                        }
                    }
                    else {
                        std::cerr << "L'Hopital's Rule: Applying derivative (application #" << recursion_depth + 1 << ").\n";
                    }

                    NodePtr new_f = make_op("/", d_num, d_den);
                    return analytic_limit(new_f, var, val_point, recursion_depth + 1);
                }
            }
            else {
                std::cerr << "DEBUG: Numerator or denominator could not be fully substituted for L'Hopital's check. Skipping L'Hopital's rule.\n";
            }
        }

        // 4. Special handling for limits at infinity for rational functions (leading coefficients)
        if (val_point.type == LimitPointType::PositiveInfinity || val_point.type == LimitPointType::NegativeInfinity) {
            if (f->type == FuncType::Div) {
                NodePtr numerator = f->children[0];
                NodePtr denominator = f->children[1];

                double num_coeff = 0.0;
                double den_coeff = 0.0;

                std::optional<double> num_deg_opt = symbolic_polynomial_analysis::get_polynomial_degree_and_coeff(numerator, var, num_coeff);
                std::optional<double> den_deg_opt = symbolic_polynomial_analysis::get_polynomial_degree_and_coeff(denominator, var, den_coeff);

                if (num_deg_opt && den_deg_opt) {
                    double num_deg = *num_deg_opt;
                    double den_deg = *den_deg_opt;

                    std::cerr << "DEBUG: Rational function at infinity: num_deg=" << num_deg << ", den_deg=" << den_deg
                        << ", num_coeff=" << num_coeff << ", den_coeff=" << den_coeff << "\n";

                    if (num_deg > den_deg) {
                        if (std::fabs(den_coeff) > EPSILON_LIMIT) {
                            if (val_point.type == LimitPointType::NegativeInfinity && std::fmod(num_deg - den_deg, 2.0) != 0.0) {
                                return make_const(std::copysign(infinity, -(num_coeff / den_coeff)));
                            }
                            return make_const(std::copysign(infinity, num_coeff / den_coeff));
                        }
                        else {
                            std::cerr << "DEBUG: Denominator leading coefficient is zero for rational function at infinity. Falling back.\n";
                        }
                    }
                    else if (num_deg < den_deg) {
                        return make_const(0.0);
                    }
                    else { // Degrees are equal
                        if (std::fabs(den_coeff) > EPSILON_LIMIT) {
                            return make_const(num_coeff / den_coeff);
                        }
                        else {
                            std::cerr << "DEBUG: Denominator leading coefficient is zero for rational function at infinity. Falling back.\n";
                        }
                    }
                }
                else {
                    std::cerr << "DEBUG: Could not determine polynomial degrees for rational function limit at infinity.\n";
                }
            }
        }

        // 5. If all fails, return nullopt.
        std::cerr << "DEBUG: All analytic methods failed for current expression. Returning nullopt.\n";
        return std::nullopt;
    }
} // namespace limit_analytic