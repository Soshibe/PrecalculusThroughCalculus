import symbolic_function_parser;  // your parser module
import limit_analytic;            // your analytic limits module

#include <iostream>
#include <string>
#include <optional>
#include <stdexcept> // For catching exceptions
#include <cmath>     // For std::isinf, std::isnan

// Assuming NodePtr and FuncType are accessible via 'import symbolic_function_parser;'
// (which implicitly brings in NodePtr, FuncType, make_const, funcTypeToString etc., if exported)

// to_string function - UPDATED to use funcTypeToString for clearer output
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
        // For binary operations, use the funcTypeToString for the operator symbol
        return "(" + to_string(node->children[0]) + " " + limit_analytic::funcTypeToString(node->type) + " " + to_string(node->children[1]) + ")";
    case FuncType::Sin:
    case FuncType::Cos:
    case FuncType::Tan:
    case FuncType::Log:
    case FuncType::Ln:
    case FuncType::Abs:
    case FuncType::Exp:
    case FuncType::AddSubFunc: // Assuming these are unary functions
    case FuncType::SubAddFunc: // Assuming these are unary functions
        if (!node->children.empty()) {
            return limit_analytic::funcTypeToString(node->type) + "(" + to_string(node->children[0]) + ")";
        }
        return limit_analytic::funcTypeToString(node->type) + "()"; // Fallback for no args
    default: return "UNKNOWN_NODE_TYPE"; // Should not happen with exhaustive enum
    }
}


int main() {
    std::cout << "Calculus 1.4 Limit Evaluator\n";

    struct TestCase {
        std::string expr;
        std::string var;
        std::string limit_point_str; // String representation of the limit point
    };

    TestCase tests[] = {
        { "x^2 + 3*x + 2", "x", "2" },               // Direct substitution: 4 + 6 + 2 = 12
        { "(x^2 - 4)/(x-2)", "x", "2" },             // L'Hopital's: (2x)/1 -> 4
        { "sin(x)/x", "x", "0" },                    // L'Hopital's: cos(x)/1 -> 1
        { "(1 - cos(x))/x^2", "x", "0" },           // L'Hopital's (twice): (sin(x))/(2x) -> cos(x)/2 -> 1/2
        { "1/x", "x", "0+" },                        // One-sided limit, handled by copysign in try_substitute
        { "1/x", "x", "0-" },                        // One-sided limit, handled by copysign in try_substitute
        { "x", "x", "inf" },                         // +inf
        { "1/x", "x", "inf" },                       // 0
        { "x^2", "x", "-inf" },                      // +inf
        { "(exp(x) - 1)/x", "x", "0" },              // L'Hopital's: exp(x)/1 -> 1
        { "(ln(x))/(x-1)", "x", "1" },               // L'Hopital's: (1/x)/1 -> 1
        { "(x^3 - 8)/(x-2)", "x", "2" },             // L'Hopital's: (3x^2)/1 -> 12
        { "x/(x+1)", "x", "inf" },                   // Leading terms: x/x -> 1
        { "exp(x)/x^2", "x", "inf" },                // L'Hopital's (twice): exp(x)/(2x) -> exp(x)/2 -> inf
        { "ln(x)/x", "x", "inf" },                   // L'Hopital's: (1/x)/1 -> 0
        { "(x^2 + 2*x + 1)/(x^2 - 3*x + 2)", "x", "inf" } // Ratio of leading coeffs: 1
    };

    for (const auto& test : tests) {
        std::cout << "\n--- Evaluating limit of: " << test.expr
            << " as " << test.var << " -> " << test.limit_point_str << " ---\n";
        try {
            // Parse the expression string into a NodePtr
            NodePtr root = parse_function_string(test.expr);
            if (!root) {
                std::cout << "Parse error for expression: " << test.expr << "\n";
                continue;
            }

            // Call the analytic_limit function (using the overload that takes string for limit point)
            // This overload internally calls limit_analytic::parse_limit_point_string
            auto limit_opt = limit_analytic::analytic_limit(root, test.var, limit_analytic::parse_limit_point_string(test.limit_point_str));
            

            if (limit_opt) {
                auto limit_node = *limit_opt;
                if (limit_node->type == FuncType::Constant) {
                    // Check for infinity values
                    if (std::isinf(limit_node->value)) {
                        if (limit_node->value > 0) {
                            std::cout << "Limit = +Infinity\n";
                        }
                        else {
                            std::cout << "Limit = -Infinity\n";
                        }
                    }
                    else if (std::isnan(limit_node->value)) {
                        std::cout << "Limit = Undefined (NaN)\n";
                    }
                    else {
                        std::cout << "Limit = " << limit_node->value << "\n";
                    }
                }
                else {
                    // This case indicates that the limit function returned a symbolic expression
                    // that is not a constant, meaning it couldn't fully resolve to a number or infinity.
                    std::cout << "Limit evaluated to a symbolic expression (not a constant): " << to_string(limit_node) << "\n";
                }
            }
            else {
                std::cout << "Limit could not be determined analytically (returned nullopt).\n";
            }
        }
        catch (const std::exception& ex) {
            std::cout << "Exception: " << ex.what() << "\n";
        }
        std::cout << "----------------------\n";
    }

    std::cin.get();
    return 0;
}