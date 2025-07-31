// symbolic_polynomial_analysis.cpp
export module symbolic_polynomial_analysis; // Associates this compilation unit with the module

import symbolic_function_parser; // Accesses NodePtr, FuncType, make_const, simplify_recursive, are_equal, flatten_and_extract_terms
import <optional>;
import <string>;
import <cmath>; // For std::fabs, std::round, std::pow
import <map>;   // For std::map in flatten_and_extract_terms
import <vector>; // For std::vector
import <limits>; // For std::numeric_limits

namespace symbolic_polynomial_analysis {

    // A small epsilon for floating point comparisons, consistent with symbolic_function_parser's EPSILON
    constexpr double EPSILON_POLY = 1e-9; // Using the same value as symbolic_function_parser::EPSILON

    // NodePtr comparator for std::map keys (strict weak ordering)
    // This is needed if `flatten_and_extract_terms` uses std::map
    // It should be consistent with symbolic_function_parser's are_equal for NodePtr comparison.
    // For simplicity, we'll assume symbolic_function_parser::are_equal is sufficient for this comparison.
    // If not, a dedicated comparator might be needed here or imported.
    struct NodePtrCmpForPolyAnalysis {
        bool operator()(NodePtr const& a, NodePtr const& b) const {
            if (!a && !b) return false; // Equal
            if (!a) return true; // nullptr < non-nullptr
            if (!b) return false; // non-nullptr > nullptr

            // Compare by type first
            if (a->type != b->type) {
                return static_cast<int>(a->type) < static_cast<int>(b->type);
            }

            // Then by value for constants
            if (a->type == FuncType::Constant) {
                // Return true if 'a' is numerically less than 'b'
                return a->value < b->value - EPSILON_POLY;
            }

            // Then by name for variables and functions (if names differ)
            if (a->name != b->name) {
                return a->name < b->name;
            }

            // Then by number of children
            if (a->children.size() != b->children.size()) {
                return a->children.size() < b->children.size();
            }

            // Then recursively by children (lexicographical comparison)
            for (size_t i = 0; i < a->children.size(); ++i) {
                if ((*this)(a->children[i], b->children[i])) return true; // a[i] < b[i]
                if ((*this)(b->children[i], a->children[i])) return false; // b[i] < a[i] (means a[i] is not < b[i])
            }
            return false; // Nodes are equivalent
        }
    };

    export std::optional<double> get_polynomial_degree_and_coeff(NodePtr node, const std::string& var_name, double& out_coeff) {
        if (!node) {
            out_coeff = 0.0;
            return 0.0; // A null node can be considered 0, degree 0
        }

        // It's crucial that `simplify` is effective for `collect_like_terms`.
        // Ensure the node is simplified for accurate analysis.
        // We call symbolic_function_parser::simplify_recursive as it's the exported simplify function.
        node = simplify_recursive(node);

        // Case 1: Constant
        if (node->type == FuncType::Constant) {
            out_coeff = node->value;
            return 0.0; // Degree of a constant is 0
        }

        // Case 2: Variable (e.g., 'x', 'y')
        if (node->type == FuncType::Variable) {
            if (node->name == var_name) {
                out_coeff = 1.0;
                return 1.0; // Degree of 'x' is 1
            }
            else {
                out_coeff = 0.0; // Treat other variables as constants
                return 0.0; // Degree of 'y' (if var_name is 'x') is 0
            }
        }

        // Case 3: Negation (e.g., -expr)
        if (node->type == FuncType::Neg) {
            double inner_coeff = 0.0;
            std::optional<double> inner_degree = get_polynomial_degree_and_coeff(node->children[0], var_name, inner_coeff);
            if (inner_degree.has_value()) {
                out_coeff = -inner_coeff;
                return inner_degree;
            }
            return std::nullopt; // Not a polynomial if child isn't
        }

        // Case 4: Multiplication (e.g., C * V, V * C, V1 * V2, etc.)
        if (node->type == FuncType::Mul) {
            NodePtr left = node->children[0];
            NodePtr right = node->children[1];

            double coeff_l = 0.0, coeff_r = 0.0;
            std::optional<double> deg_l = get_polynomial_degree_and_coeff(left, var_name, coeff_l);
            std::optional<double> deg_r = get_polynomial_degree_and_coeff(right, var_name, coeff_r);

            if (deg_l.has_value() && deg_r.has_value()) {
                // Product of two polynomials is a polynomial
                out_coeff = coeff_l * coeff_r;
                return deg_l.value() + deg_r.value();
            }
            // If one is a constant and the other is not a polynomial in `var_name`
            // (e.g., 2 * sin(x) where var_name is x)
            // It's only a polynomial if the non-polynomial part evaluates to a constant with respect to var_name.
            // Let's re-evaluate more strictly: if *any* part is not a polynomial, the product is not.
            return std::nullopt; // Not a polynomial if any part is not a polynomial
        }

        // Case 5: Division (e.g., expr / expr)
        // A division is generally a polynomial only if the denominator is a non-zero constant.
        if (node->type == FuncType::Div) {
            NodePtr numerator = node->children[0];
            NodePtr denominator = node->children[1];

            double den_coeff = 0.0;
            std::optional<double> den_deg = get_polynomial_degree_and_coeff(denominator, var_name, den_coeff);

            // Denominator must be a constant (degree 0) and non-zero
            if (den_deg.has_value() && std::fabs(den_deg.value()) < EPSILON_POLY && std::fabs(den_coeff) > EPSILON_POLY) {
                double num_coeff = 0.0;
                std::optional<double> num_deg = get_polynomial_degree_and_coeff(numerator, var_name, num_coeff);
                if (num_deg.has_value()) {
                    out_coeff = num_coeff / den_coeff;
                    return num_deg;
                }
            }
            return std::nullopt; // Not a polynomial
        }

        // Case 6: Power (e.g., expr ^ N)
        // A power is a polynomial if the base is `var_name` and exponent is a non-negative integer constant,
        // or if the base is a constant and exponent is a non-negative integer constant (resulting in a constant),
        // or if the base is a polynomial and exponent is a non-negative integer constant.
        if (node->type == FuncType::Pow) {
            NodePtr base = node->children[0];
            NodePtr exponent = node->children[1];

            if (exponent->type == FuncType::Constant) {
                double exp_val = exponent->value;
                // Exponent must be a non-negative integer for it to be a polynomial term
                if (exp_val >= 0 && std::fabs(exp_val - std::round(exp_val)) < EPSILON_POLY) {
                    double base_coeff = 0.0;
                    std::optional<double> base_deg = get_polynomial_degree_and_coeff(base, var_name, base_coeff);
                    if (base_deg.has_value()) {
                        // For (C*x^n)^m, the new coefficient is C^m and new degree is n*m
                        // This handles cases like (2*x)^3 = 8*x^3
                        out_coeff = std::pow(base_coeff, static_cast<int>(std::round(exp_val)));
                        return base_deg.value() * exp_val;
                    }
                }
            }
            return std::nullopt; // Not a polynomial power
        }

        // Case 7: Addition/Subtraction (e.g., expr1 + expr2)
        // The sum/difference of polynomials is a polynomial. The degree is the max of the children's degrees.
        // The coefficient is trickier: it's the coefficient of the highest degree term after collection.
        if (node->type == FuncType::Add || node->type == FuncType::Sub) {
            // We rely on `simplify_recursive` having `collect_like_terms` working well.
            // If terms are collected, an expression like `2x + 3x` would already be `5x`.
            // Let's re-use `flatten_and_extract_terms` to find the highest degree from the simplified form.

            std::vector<std::pair<double, NodePtr>> terms = flatten_and_extract_terms(node);
            double max_degree = -1.0; // Use -1 to represent "not a polynomial" or "no terms found"

            // Map to hold degrees and their accumulated coefficients
            // Use NodePtrCmpForPolyAnalysis to correctly group terms by their base (e.g., x^2, x, 1)
            std::map<double, double> degree_to_coeff_sum;

            for (const auto& term_pair : terms) {
                double term_coeff_from_flatten = term_pair.first;
                NodePtr base_node = term_pair.second;

                double current_term_base_coeff_val = 0.0;
                std::optional<double> current_term_degree_val = std::nullopt;

                // Analyze the base_node (e.g., x^2, x, 1, y^3 if var_name is x)
                if (base_node->type == FuncType::Variable && base_node->name == var_name) {
                    current_term_degree_val = 1.0;
                    current_term_base_coeff_val = 1.0;
                }
                else if (base_node->type == FuncType::Constant) {
                    current_term_degree_val = 0.0;
                    current_term_base_coeff_val = base_node->value;
                }
                else if (base_node->type == FuncType::Pow && base_node->children[0]->type == FuncType::Variable && base_node->children[0]->name == var_name) {
                    // Check if the exponent is a non-negative integer constant
                    NodePtr exp_node = base_node->children[1];
                    if (exp_node->type == FuncType::Constant) {
                        double exp_val = exp_node->value;
                        if (exp_val >= 0 && std::fabs(exp_val - std::round(exp_val)) < EPSILON_POLY) {
                            current_term_degree_val = exp_val;
                            current_term_base_coeff_val = 1.0; // The base itself is 'x', so its coefficient is 1
                        }
                        else {
                            out_coeff = 0.0; // Not a polynomial due to non-integer or negative exponent
                            return std::nullopt;
                        }
                    }
                    else {
                        out_coeff = 0.0; // Not a polynomial because exponent is not constant
                        return std::nullopt;
                    }
                }
                else {
                    // If the base node is not x, a constant, or x^N, it must be a constant with respect to `var_name`
                    // (e.g., 'y' when var_name is 'x', or 'sin(z)' if z is not var_name).
                    // If it contains `var_name` in a non-polynomial way (e.g., sin(x)), it's not a polynomial.
                    double temp_coeff_for_base_check = 0.0;
                    std::optional<double> base_node_deg = get_polynomial_degree_and_coeff(base_node, var_name, temp_coeff_for_base_check);

                    if (base_node_deg.has_value() && std::fabs(base_node_deg.value()) < EPSILON_POLY) {
                        // It's a constant w.r.t. var_name (degree 0)
                        current_term_degree_val = 0.0;
                        current_term_base_coeff_val = temp_coeff_for_base_check; // The 'constant' value of this base node
                    }
                    else {
                        out_coeff = 0.0; // Not a polynomial (e.g., sin(x) term, or x*y where y is not var_name)
                        return std::nullopt;
                    }
                }

                if (current_term_degree_val.has_value()) {
                    // Accumulate coefficients for each degree
                    degree_to_coeff_sum[current_term_degree_val.value()] += term_coeff_from_flatten * current_term_base_coeff_val;
                }
                else {
                    out_coeff = 0.0; // Not a polynomial
                    return std::nullopt;
                }
            }

            double final_leading_coeff = 0.0;

            // Iterate in reverse to find the highest degree (map sorts by key ascending)
            for (auto it = degree_to_coeff_sum.rbegin(); it != degree_to_coeff_sum.rend(); ++it) {
                double degree = it->first;
                double coeff_sum = it->second;

                if (std::fabs(coeff_sum) > EPSILON_POLY) { // Only consider terms with non-zero coefficients
                    max_degree = degree;
                    final_leading_coeff = coeff_sum;
                    break; // Found the highest degree non-zero term
                }
            }

            // If all coefficients summed to zero (e.g., x - x), the result is a polynomial of degree 0 with coeff 0
            if (std::fabs(final_leading_coeff) < EPSILON_POLY && max_degree == -1.0) {
                out_coeff = 0.0;
                return 0.0; // All terms cancelled out, effectively 0 constant
            }
            else if (max_degree == -1.0) { // No terms found or all terms had zero coefficients
                out_coeff = 0.0;
                return 0.0; // This case might occur if input was just '0'
            }

            out_coeff = final_leading_coeff;
            return max_degree;
        }

        // Default: If none of the above, it's not a polynomial with respect to `var_name`
        // e.g., sin(x), log(x), y^x, etc.
        out_coeff = 0.0;
        return std::nullopt;
    }

} // namespace symbolic_polynomial_analysis