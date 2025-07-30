// main.cpp
import symbolic_function_parser;
import symbolic_function_differentiation;
import limit_analytic;

#include <iostream>
import <string>;
import <vector>;
import <memory>; // For NodePtr
import <optional>;

// Function to print the symbolic tree (useful for debugging)
std::string print_node(NodePtr node) {
    if (!node) return "nullptr";

    std::string result = "";
    switch (node->type) {
    case FuncType::Constant:
        result += std::to_string(node->value);
        break;
    case FuncType::Variable:
        result += node->name;
        break;
    case FuncType::Add:
    case FuncType::Sub:
    case FuncType::Mul:
    case FuncType::Div:
    case FuncType::Pow:
        result += "(" + print_node(node->children[0]) + node->name + print_node(node->children[1]) + ")";
        break;
    case FuncType::Neg:
        result += "(-" + print_node(node->children[0]) + ")";
        break;
    case FuncType::Sin: case FuncType::Cos: case FuncType::Tan:
    case FuncType::Ln: case FuncType::Log: case FuncType::Abs:
    case FuncType::Exp:
        result += node->name + "(" + print_node(node->children[0]) + ")";
        break;
    case FuncType::AddSubFunc:
        result += "addsub(" + print_node(node->children[0]) + ")";
        break;
    case FuncType::SubAddFunc:
        result += "subadd(" + print_node(node->children[0]) + ")";
        break;
    case FuncType::UnknownFunc:
        result += "UNKNOWN(" + node->name + ")";
        if (!node->children.empty()) {
            result += "(";
            for (size_t i = 0; i < node->children.size(); ++i) {
                result += print_node(node->children[i]);
                if (i < node->children.size() - 1) result += ",";
            }
            result += ")";
        }
        break;
    }
    return result;
}

int main() {
    std::cout << "Symbolic Expression Calculator Example" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    // --- Example 1: Basic Parsing and Simplification ---
    std::string expr1_str = "2*x + 3*x - 5";
    std::cout << "\nExpression 1: " << expr1_str << std::endl;
    try {
        NodePtr expr1_parsed = parse_function_string(expr1_str);
        std::cout << "Parsed (raw): " << print_node(expr1_parsed) << std::endl;
        NodePtr expr1_simplified = simplify(expr1_parsed);
        std::cout << "Simplified: " << print_node(expr1_simplified) << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // --- Example 2: Derivative ---
    std::string expr2_str = "sin(x^2) + 3*x";
    std::string var2 = "x";
    std::cout << "\nExpression 2: " << expr2_str << std::endl;
    std::cout << "Differentiate with respect to: " << var2 << std::endl;
    try {
        NodePtr expr2_parsed = parse_function_string(expr2_str);
        NodePtr d_expr2 = derivative(expr2_parsed, var2);
        std::cout << "Derivative (raw): " << print_node(d_expr2) << std::endl;
        NodePtr d_expr2_simplified =  simplify(d_expr2);
        std::cout << "Derivative (simplified): " << print_node(d_expr2_simplified) << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // --- Example 3: Limit using L'Hôpital's Rule (0/0 form) ---
    // Limit as x->0 of sin(x)/x
    std::string expr3_str = "sin(x)/x";
    std::string var3 = "x";
    double val3 = 0.0;
    std::cout << "\nLimit Example 1: limit as " << var3 << "->" << val3 << " of " << expr3_str << std::endl;
    try {
        NodePtr expr3_parsed =  parse_function_string(expr3_str);
        std::optional<NodePtr> limit3_result = analytic_limit(expr3_parsed, var3, val3);
        if (limit3_result.has_value()) {
            std::cout << "Limit result: " << print_node(limit3_result.value()) << std::endl;
        }
        else {
            std::cout << "Could not evaluate limit analytically." << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // --- Example 4: Limit with direct substitution ---
    // Limit as x->2 of x^2 + 1
    std::string expr4_str = "x^2 + 1";
    std::string var4 = "x";
    double val4 = 2.0;
    std::cout << "\nLimit Example 2: limit as " << var4 << "->" << val4 << " of " << expr4_str << std::endl;
    try {
        NodePtr expr4_parsed =  parse_function_string(expr4_str);
        std::optional<NodePtr> limit4_result = analytic_limit(expr4_parsed, var4, val4);
        if (limit4_result.has_value()) {
            std::cout << "Limit result: " << print_node(limit4_result.value()) << std::endl;
        }
        else {
            std::cout << "Could not evaluate limit analytically." << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // --- Example 5: Limit with AddSubFunc (should simplify to 0 if argument is 0) ---
    // Limit as x->0 of addsub(x)
    std::string expr5_str = "addsub(x)";
    std::string var5 = "x";
    double val5 = 0.0;
    std::cout << "\nLimit Example 3: limit as " << var5 << "->" << val5 << " of " << expr5_str << std::endl;
    try {
        NodePtr expr5_parsed =  parse_function_string(expr5_str);
        std::optional<NodePtr> limit5_result =  analytic_limit(expr5_parsed, var5, val5);
        if (limit5_result.has_value()) {
            std::cout << "Limit result: " << print_node(limit5_result.value()) << std::endl;
        }
        else {
            std::cout << "Could not evaluate limit analytically." << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // --- Example 6: Limit that might lead to a non-numerical result or error (e.g., variable in limit result) ---
    // Limit as x->0 of x+y
    std::string expr6_str = "(5/x)^2+3+5x";
    std::string var6 = "x";
    double val6 = 0.0;
    std::cout << "\nLimit Example 4: limit as " << var6 << "->" << val6 << " of " << expr6_str << std::endl;
    try {
        NodePtr expr6_parsed =  parse_function_string(expr6_str);
        std::optional<NodePtr> limit6_result =  analytic_limit(expr6_parsed, var6, val6);
        if (limit6_result.has_value()) {
            std::cout << "Limit result: " << print_node(limit6_result.value()) << std::endl;
        }
        else {
            std::cout << "Could not evaluate limit analytically (expected for symbolic result)." << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    std::cin.get();
    return 0;
}