// symbolic_function_parser.ixx
export module symbolic_function_parser;

import <string>;
import <vector>;
import <memory>;
import <optional>;
import <cctype>;
import <map>;
import <stdexcept>;
import <cmath>;
import <numeric>;
import <algorithm>; // For std::sort (if needed, not directly used here but good practice)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

export enum class FuncType {
    Constant,
    Variable,
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Neg,
    Sin,
    Cos,
    Tan,
    Log,
    Ln,
    Abs,
    Exp,
    AddSubFunc,
    SubAddFunc,
    UnknownFunc
};

export struct SymbolicNode {
    FuncType type;
    std::string name;
    double value = 0.0;
    std::vector<std::shared_ptr<SymbolicNode>> children;
    bool is_primary = false;
    std::string original_expression; // Consider if this is truly needed or can be removed
};

export using NodePtr = std::shared_ptr<SymbolicNode>;

export NodePtr make_op(const std::string& op, NodePtr lhs, NodePtr rhs) {
    FuncType type = FuncType::UnknownFunc;
    if (op == "+") type = FuncType::Add;
    else if (op == "-") type = FuncType::Sub;
    else if (op == "*") type = FuncType::Mul;
    else if (op == "/") type = FuncType::Div;
    else if (op == "^") type = FuncType::Pow;
    else throw std::runtime_error("Unknown binary operator: " + op);
    return std::make_shared<SymbolicNode>(SymbolicNode{ type, op, 0.0, {std::move(lhs), std::move(rhs)} }); // Store op as name for debugging
}

export NodePtr make_func(FuncType ft, const std::string& name, std::vector<NodePtr> args) {
    return std::make_shared<SymbolicNode>(SymbolicNode{ ft, name, 0.0, std::move(args) });
}

export NodePtr make_var(const std::string& name) {
    return std::make_shared<SymbolicNode>(SymbolicNode{ FuncType::Variable, name });
}

export NodePtr make_const(double val) {
    return std::make_shared<SymbolicNode>(SymbolicNode{ FuncType::Constant, "", val });
}

export NodePtr make_unary(FuncType ft, NodePtr arg) {
    std::string name = ""; // Unary ops usually don't have a 'name' like sin/cos
    if (ft == FuncType::Neg) name = "-"; // For -x
    return std::make_shared<SymbolicNode>(SymbolicNode{ ft, name, 0.0, {std::move(arg)} });
}



    static constexpr double EPSILON = 1e-9; // Increased precision for double comparisons

    const std::map<std::string, FuncType> func_lookup = {
        {"sin", FuncType::Sin}, {"cos", FuncType::Cos}, {"tan", FuncType::Tan},
        {"log", FuncType::Log}, {"ln", FuncType::Ln}, {"abs", FuncType::Abs},
        {"exp", FuncType::Exp},
        {"addsub", FuncType::AddSubFunc},
        {"subadd", FuncType::SubAddFunc}
    };

    struct Token {
        enum class Type { Number, Identifier, Operator, LParen, RParen, End } type;
        std::string text;
        std::string to_string() const {
            switch (type) {
            case Type::Number: return "Number(" + text + ")";
            case Type::Identifier: return "Identifier(\"" + text + "\")";
            case Type::Operator: return "Operator('" + text + "')";
            case Type::LParen: return "LParen";
            case Type::RParen: return "RParen";
            case Type::End: return "End";
            }
            return "Unknown";
        }
    };

    class Tokenizer {
        std::string_view expr;
        size_t pos = 0;
    public:
        Tokenizer(std::string_view s) : expr(s) {}
        Token next_token() {
            while (pos < expr.size() && std::isspace(static_cast<unsigned char>(expr[pos]))) ++pos;
            if (pos >= expr.size()) return { Token::Type::End, "" };

            char c = expr[pos];
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
                size_t start = pos;
                while (pos < expr.size() && (std::isdigit(static_cast<unsigned char>(expr[pos])) || expr[pos] == '.')) ++pos;
                return { Token::Type::Number, std::string(expr.substr(start, pos - start)) };
            }
            if (std::isalpha(static_cast<unsigned char>(c))) {
                size_t start = pos;
                while (pos < expr.size() && (std::isalnum(static_cast<unsigned char>(expr[pos])) || expr[pos] == '_')) ++pos;
                return { Token::Type::Identifier, std::string(expr.substr(start, pos - start)) };
            }
            if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
                ++pos;
                return { Token::Type::Operator, std::string(1, c) };
            }
            if (c == '(') { ++pos; return { Token::Type::LParen, "(" }; }
            if (c == ')') { ++pos; return { Token::Type::RParen, ")" }; }

            throw std::runtime_error("Invalid character in expression: '" + std::string(1, c) + "' at position " + std::to_string(pos));
        }
        Token peek_next_token() const {
            Tokenizer temp_tokenizer(expr.substr(pos));
            return temp_tokenizer.next_token();
        }
    };

    class Parser {
        Tokenizer tokenizer;
        Token current_token;
        void advance() { current_token = tokenizer.next_token(); }
        [[nodiscard]] bool match(Token::Type type, const std::string& text = "") {
            return current_token.type == type && (text.empty() || current_token.text == text);
        }
        NodePtr parse_expression(int min_prec = 0) {
            NodePtr lhs = parse_factor();

            while (true) {
                int prec = get_precedence(current_token);
                if (prec < min_prec) break;

                std::string op = current_token.text;
                advance(); // Consume the operator

                NodePtr rhs = parse_expression(prec + (op == "^" ? 0 : 1)); // Right-associativity for power
                lhs = make_op(op, std::move(lhs), std::move(rhs));
            }
            return lhs;
        }

        NodePtr parse_factor() {
            NodePtr node = parse_atomic_expression(); // Parses '2', current_token is now '*'

            // After parsing an atomic expression, current_token points to the *next* token.
            // If that token is an operator, RParen, or End, it's not implicit multiplication.
            // It should be handled by parse_expression's loop.
            while (current_token.type != Token::Type::Operator &&
                current_token.type != Token::Type::RParen &&
                current_token.type != Token::Type::End)
            {
                // If we reach here, current_token is NOT an operator, RParen, or End.
                // It must be a Number, Identifier, or LParen, which signifies implicit multiplication.
                // Example: '2x', 'x(y+z)', 'sin(x)y'
                NodePtr rhs_atomic = parse_atomic_expression(); // This will consume current_token and advance
                node = make_op("*", std::move(node), std::move(rhs_atomic));
            }
            return node;
        }

        NodePtr parse_atomic_expression() {
            if (match(Token::Type::Operator, "-")) {
                advance();
                return make_unary(FuncType::Neg, parse_atomic_expression());
            }
            if (match(Token::Type::Number)) {
                double val = std::stod(current_token.text);
                advance();
                return make_const(val);
            }
            if (match(Token::Type::Identifier)) {
                std::string id = current_token.text;
                advance();
                if (id == "pi" || id == "PI") {
                    return make_const(M_PI);
                }
                if (match(Token::Type::LParen)) {
                    advance();
                    // Handle comma-separated arguments for future multi-arg functions
                    // For now, only single argument expected.
                    auto arg = parse_expression();
                    if (!match(Token::Type::RParen)) {
                        throw std::runtime_error("Expected ')' after function argument, got " + current_token.to_string());
                    }
                    advance();
                    auto it = func_lookup.find(id);
                    FuncType ft = (it != func_lookup.end()) ? it->second : FuncType::UnknownFunc;
                    return make_func(ft, id, { std::move(arg) });
                }
                return make_var(id);
            }
            if (match(Token::Type::LParen)) {
                advance();
                auto node = parse_expression();
                if (!match(Token::Type::RParen)) {
                    throw std::runtime_error("Expected ')' after expression in parentheses, got " + current_token.to_string());
                }
                advance();
                return node;
            }
            throw std::runtime_error("Unexpected token while parsing atomic expression: " + current_token.to_string());
        }
        static int get_precedence(const Token& tok) {
            if (tok.type != Token::Type::Operator) return -1;
            if (tok.text == "+" || tok.text == "-") return 1;
            if (tok.text == "*" || tok.text == "/") return 2;
            if (tok.text == "^") return 3;
            return -1;
        }
    public:
        explicit Parser(std::string_view input) : tokenizer(input) {
            advance(); // Initialize current_token with the first token
        }
        NodePtr parse() {
            NodePtr result = parse_expression();
            if (current_token.type != Token::Type::End) {
                throw std::runtime_error("Unexpected input after parsing: " + current_token.text);
            }
            return result;
        }
    };

    // Forward declarations for mutual recursion in are_equal and simplify_recursive
    export bool are_equal(NodePtr n1, NodePtr n2);
    NodePtr simplify_recursive(NodePtr node);

    // Helper for deep comparison of NodePtrs
    bool are_equal_vectors(const std::vector<NodePtr>& v1, const std::vector<NodePtr>& v2) {
        if (v1.size() != v2.size()) return false;
        for (size_t i = 0; i < v1.size(); ++i) {
            // Check if both elements are non-null before dereferencing
            if (!v1[i] || !v2[i]) {
                if (v1[i] != v2[i]) return false; // One is null, other isn't
                continue; // Both are null, continue to next element
            }
            if (!are_equal(v1[i], v2[i])) { // Recursive call
                return false;
            }
        }
        return true;
    }

    bool are_equal(NodePtr n1, NodePtr n2) {
        if (!n1 && !n2) return true;
        if (!n1 || !n2) return false;

        if (n1->type != n2->type) return false;

        switch (n1->type) {
        case FuncType::Constant:
            return std::fabs(n1->value - n2->value) < EPSILON;
        case FuncType::Variable:
            return n1->name == n2->name;
            // All other types rely on name (if applicable) and children equality
        case FuncType::Add:
        case FuncType::Sub:
        case FuncType::Mul:
        case FuncType::Div:
        case FuncType::Pow:
        case FuncType::Neg:
        case FuncType::Sin:
        case FuncType::Cos:
        case FuncType::Tan:
        case FuncType::Log:
        case FuncType::Ln:
        case FuncType::Abs:
        case FuncType::Exp:
        case FuncType::AddSubFunc:
        case FuncType::SubAddFunc:
        case FuncType::UnknownFunc:
            // Compare names for operators like "+", "-", and named functions
            if (n1->name != n2->name) return false;
            return are_equal_vectors(n1->children, n2->children);
        default:
            return false;
        }
    }

    // NodePtr comparator for std::map keys (strict weak ordering)
    struct NodePtrCmp {
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
                return a->value < b->value - EPSILON;
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


    bool is_constant_val(NodePtr node, double val) {
        return node && node->type == FuncType::Constant && std::fabs(node->value - val) < EPSILON;
    }

    // Helper: Flattens an expression into a list of (coefficient, base) pairs.
    // Handles nesting of Add/Sub for a "sum of terms" view.
    // The `base_node` is the non-constant, non-negative part of a term.
    // For `5`, it's `(5.0, make_const(1.0))`
    // For `2*x`, it's `(2.0, x_NodePtr)`
    // For `x`, it's `(1.0, x_NodePtr)`
    // For `-x`, it's `(-1.0, x_NodePtr)`
    // For `-(2*x)`, it's `(-2.0, x_NodePtr)`
    std::vector<std::pair<double, NodePtr>> flatten_and_extract_terms(NodePtr node) {
        std::vector<std::pair<double, NodePtr>> terms;
        if (!node) return terms;

        if (node->type == FuncType::Constant) {
            terms.push_back({ node->value, make_const(1.0) }); // Constant terms have value as coeff, base is 1
            return terms;
        }
        if (node->type == FuncType::Neg) { // Handle -A as -1 * A, then flatten A
            auto inner_terms = flatten_and_extract_terms(node->children[0]);
            for (auto& term_pair : inner_terms) {
                term_pair.first *= -1.0;
            }
            return inner_terms;
        }
        if (node->type == FuncType::Mul) { // Handle C*X or X*C terms
            NodePtr left = node->children[0];
            NodePtr right = node->children[1];
            if (left->type == FuncType::Constant) {
                // If it's a constant multiplied by something (e.g., 2*x)
                auto sub_terms = flatten_and_extract_terms(right); // Flatten the 'something'
                for (auto& term_pair : sub_terms) {
                    term_pair.first *= left->value; // Multiply its coefficient by our constant
                }
                return sub_terms;
            }
            if (right->type == FuncType::Constant) {
                // If it's something multiplied by a constant (e.g., x*2)
                auto sub_terms = flatten_and_extract_terms(left); // Flatten the 'something'
                for (auto& term_pair : sub_terms) {
                    term_pair.first *= right->value; // Multiply its coefficient by our constant
                }
                return sub_terms;
            }
        }
        // For Add and Sub, recursively flatten children.
        if (node->type == FuncType::Add) {
            auto left_terms = flatten_and_extract_terms(node->children[0]);
            auto right_terms = flatten_and_extract_terms(node->children[1]);
            terms.insert(terms.end(), left_terms.begin(), left_terms.end());
            terms.insert(terms.end(), right_terms.begin(), right_terms.end());
            return terms;
        }
        if (node->type == FuncType::Sub) { // A - B => A + (-B) => A + (-1)*B
            auto left_terms = flatten_and_extract_terms(node->children[0]);
            auto right_terms_negated = flatten_and_extract_terms(node->children[1]);
            for (auto& term_pair : right_terms_negated) {
                term_pair.first *= -1.0;
            }
            terms.insert(terms.end(), left_terms.begin(), left_terms.end());
            terms.insert(terms.end(), right_terms_negated.begin(), right_terms_negated.end());
            return terms;
        }

        // Default: it's a base term with a coefficient of 1
        terms.push_back({ 1.0, node });
        return terms;
    }

    // This function assumes its input `node` has its children already recursively simplified.
    // It groups terms and reconstructs the expression.
    NodePtr collect_like_terms(NodePtr node, bool& changed) {
        if (!node || (node->type != FuncType::Add && node->type != FuncType::Sub)) {
            // This function should only be called on Add or Sub nodes.
            return node;
        }

        std::map<NodePtr, double, NodePtrCmp> grouped_terms;
        bool local_changed = false; // Tracks changes within this specific call

        // Flatten the current expression into its constituent terms (coeff, base)
        std::vector<std::pair<double, NodePtr>> flat_terms = flatten_and_extract_terms(node);

        for (const auto& term_pair : flat_terms) {
            double coeff = term_pair.first;
            NodePtr base = term_pair.second;

            // If the base is a constant '1.0', it means the term itself was originally a constant.
            // Sum its coefficient directly into a special constant sum.
            if (base && base->type == FuncType::Constant && std::fabs(base->value - 1.0) < EPSILON) {
                grouped_terms[make_const(1.0)] += coeff; // Use make_const(1.0) as a canonical key for all pure constants
            }
            else {
                grouped_terms[base] += coeff;
            }
            local_changed = true; // Any flattening/collection implies a potential change
        }

        NodePtr result_node = nullptr;

        // Reconstruct the expression from combined terms
        for (const auto& pair : grouped_terms) {
            NodePtr base = pair.first;
            double coeff = pair.second;

            if (std::fabs(coeff) < EPSILON) { // If coefficient became zero after combining (e.g., x - x)
                local_changed = true;
                continue;
            }

            NodePtr term_to_add;
            if (base->type == FuncType::Constant && std::fabs(base->value - 1.0) < EPSILON) {
                // This is the accumulated constant sum
                term_to_add = make_const(coeff); // The coefficient IS the constant value
            }
            else if (std::fabs(coeff - 1.0) < EPSILON) {
                term_to_add = base;
            }
            else if (std::fabs(coeff - (-1.0)) < EPSILON) {
                term_to_add = make_unary(FuncType::Neg, base);
            }
            else {
                term_to_add = make_op("*", make_const(coeff), base);
            }

            if (!result_node) {
                result_node = term_to_add;
            }
            else {
                result_node = make_op("+", result_node, term_to_add);
            }
        }

        // Handle case where all terms cancelled out, result_node might be null
        if (!result_node) {
            result_node = make_const(0.0);
            local_changed = true;
        }

        // Only set the global 'changed' flag if an actual structural change occurred
        if (local_changed && !are_equal(node, result_node)) {
            changed = true;
        }
        else {
            // If local_changed was true but no overall structural change,
            // it means the input node was already in the simplified form or
            // the changes were minor and 'are_equal' caught them.
            // In this case, return the original node if no actual change detected
            return node;
        }

        return result_node;
    }


    NodePtr simplify_once(NodePtr node, bool& changed) {
        if (!node) return nullptr;

        // Step 1: Recursively simplify all children first
        std::vector<NodePtr> simplified_children;
        simplified_children.reserve(node->children.size());
        bool children_changed = false; // Track if *any* child was changed
        for (const auto& c : node->children) {
            NodePtr simplified_c = simplify_recursive(c); // Recursive call
            if (!are_equal(simplified_c, c)) {
                children_changed = true;
            }
            simplified_children.push_back(std::move(simplified_c));
        }
        // Only update node's children if any child actually changed to avoid unnecessary re-creation
        if (children_changed) {
            node->children = std::move(simplified_children);
            changed = true; // Mark as changed because children were updated
        }

        // If the node itself is a constant or variable, it cannot be simplified further at this level,
        // unless its children were changed (which is already handled above).
        if (node->type == FuncType::Constant || node->type == FuncType::Variable) {
            return node;
        }

        // Step 2: Apply simplification rules to the current node

        // Handle specific cases for AddSubFunc and SubAddFunc:
        if (node->type == FuncType::AddSubFunc || node->type == FuncType::SubAddFunc) {
            if (node->children.size() == 1 && is_constant_val(node->children[0], 0.0)) {
                changed = true; // Rule application
                return make_const(0.0);
            }
            return node; // Return the node with potentially simplified children
        }

        // Constant folding for binary operators (Add, Sub, Mul, Div, Pow)
        if ((node->type == FuncType::Add || node->type == FuncType::Sub ||
            node->type == FuncType::Mul || node->type == FuncType::Div ||
            node->type == FuncType::Pow) &&
            node->children.size() == 2 &&
            node->children[0]->type == FuncType::Constant &&
            node->children[1]->type == FuncType::Constant)
        {
            changed = true; // This rule applies, so a change will occur
            double lhs = node->children[0]->value;
            double rhs = node->children[1]->value;
            switch (node->type) {
            case FuncType::Add: return make_const(lhs + rhs);
            case FuncType::Sub: return make_const(lhs - rhs);
            case FuncType::Mul: return make_const(lhs * rhs);
            case FuncType::Div:
                if (std::fabs(rhs) < EPSILON) throw std::runtime_error("Division by zero in simplification");
                return make_const(lhs / rhs);
            case FuncType::Pow:
                if (lhs < 0.0 && std::fabs(std::floor(rhs) - rhs) > EPSILON) // fractional exponent of negative base
                    throw std::runtime_error("Negative base with non-integer exponent in constant folding");
                return make_const(std::pow(lhs, rhs));
            default: break; // Should not happen
            }
        }

        // Unary constant folding (Neg, Sin, Cos, Tan, Ln, Log, Abs, Exp)
        if (node->children.size() == 1 && node->children[0]->type == FuncType::Constant) {
            double val = node->children[0]->value;
            switch (node->type) {
            case FuncType::Neg: changed = true; return make_const(-val);
            case FuncType::Sin: changed = true; return make_const(std::sin(val));
            case FuncType::Cos: changed = true; return make_const(std::cos(val));
            case FuncType::Tan:
                // Check for multiples of pi/2 (where tan is undefined)
                if (std::fabs(std::fmod(std::fabs(val), M_PI) - M_PI_2) < EPSILON)
                    throw std::runtime_error("tan(x) undefined at pi/2 + n*pi in constant folding");
                changed = true; return make_const(std::tan(val));
            case FuncType::Ln:
                if (val <= EPSILON) throw std::runtime_error("ln(x) undefined for x <= 0 in constant folding");
                changed = true; return make_const(std::log(val));
            case FuncType::Log:
                if (val <= EPSILON) throw std::runtime_error("log(x) undefined for x <= 0 in constant folding");
                changed = true; return make_const(std::log10(val));
            case FuncType::Abs: changed = true; return make_const(std::fabs(val));
            case FuncType::Exp: changed = true; return make_const(std::exp(val));
            default: break;
            }
        }

        // Specific algebraic simplifications (after children are simplified and constants folded)
        switch (node->type) {
        case FuncType::Add:
        case FuncType::Sub: {
            // Attempt to collect and combine like terms.
            // This function internally sets 'changed' if a transformation occurs.
            NodePtr combined_node = collect_like_terms(node, changed);
            return combined_node; // Always return the result of collect_like_terms
        }
        case FuncType::Mul: {
            NodePtr left = node->children[0];
            NodePtr right = node->children[1];
            if (is_constant_val(left, 1.0)) { changed = true; return right; }
            if (is_constant_val(right, 1.0)) { changed = true; return left; }
            if (is_constant_val(left, 0.0) || is_constant_val(right, 0.0)) {
                changed = true; return make_const(0.0);
            }
            // -1 * x = -x
            if (is_constant_val(left, -1.0)) { changed = true; return make_unary(FuncType::Neg, right); }
            if (is_constant_val(right, -1.0)) { changed = true; return make_unary(FuncType::Neg, left); }
            break;
        }
        case FuncType::Div: {
            NodePtr left = node->children[0];
            NodePtr right = node->children[1];
            if (is_constant_val(right, 1.0)) { changed = true; return left; }
            if (is_constant_val(left, 0.0) && !is_constant_val(right, 0.0)) {
                changed = true; return make_const(0.0);
            }
            if (is_constant_val(right, -1.0)) { changed = true; return make_unary(FuncType::Neg, left); }
            if (are_equal(left, right) && !is_constant_val(left, 0.0)) {
                changed = true; return make_const(1.0);
            }
            break;
        }
        case FuncType::Pow: {
            NodePtr base = node->children[0];
            NodePtr exponent = node->children[1];
            if (is_constant_val(exponent, 0.0)) { // x^0 = 1 (if x is not 0)
                if (is_constant_val(base, 0.0)) { /* 0^0 is indeterminate, not simplified to 1 */ return node; }
                changed = true; return make_const(1.0);
            }
            if (is_constant_val(exponent, 1.0)) { changed = true; return base; } // x^1 = x
            if (is_constant_val(base, 0.0) && // 0^x = 0 (if x > 0)
                exponent->type == FuncType::Constant &&
                exponent->value > EPSILON) {
                changed = true; return make_const(0.0);
            }
            if (is_constant_val(base, 1.0)) { changed = true; return make_const(1.0); } // 1^x = 1
            break;
        }
        case FuncType::Neg: {
            NodePtr child = node->children[0];
            if (child->type == FuncType::Neg) { // --x = x
                changed = true;
                return child->children[0];
            }
            break;
        }
        case FuncType::Exp: {
            NodePtr arg = node->children[0];
            if (arg->type == FuncType::Ln) { // exp(ln(x)) = x
                changed = true;
                return arg->children[0];
            }
            break;
        }
        case FuncType::Ln: {
            NodePtr arg = node->children[0];
            if (arg->type == FuncType::Exp) { // ln(exp(x)) = x
                changed = true;
                return arg->children[0];
            }
            break;
        }
        case FuncType::Sin:
        case FuncType::Cos:
        case FuncType::Tan:
        case FuncType::Log:
        case FuncType::Abs:
        case FuncType::UnknownFunc: // Handle UnknownFunc here if it has children and needs to return the updated node
        default:
            // Variables, constants, and other unknown/unhandled functions
            // are returned as is (or already handled by specific rules above).
            break;
        }

        return node; // Return the node if no rules applied at this level or it's a base type
    }

    NodePtr simplify_recursive(NodePtr node) {
        if (!node) return nullptr;
        // Base cases: If it's a constant or variable, no structural simplification possible at this level,
        // unless children were changed, which is handled in simplify_once's first step.
        // UnknownFunc also acts as a base case unless it has defined simplification rules later.
        if (node->type == FuncType::Constant || node->type == FuncType::Variable)
            return node;

        bool changed = false;
        NodePtr current = node;
        do {
            changed = false;
            NodePtr previous_state_for_comparison = current; // Capture current state for comparison

            // This is the core simplification step.
            // `current` might be replaced by a new node, or its children might be updated.
            current = simplify_once(current, changed); // Pass `current` by value, as `simplify_once` often returns a new node

            // If simplify_once didn't set 'changed' but the node structure actually changed,
            // we catch it here to ensure the loop continues.
            // This happens if a new NodePtr is returned from simplify_once without `changed` being set to true inside it.
            // `are_equal` is a deep structural comparison.
            if (!changed && !are_equal(previous_state_for_comparison, current)) {
                changed = true;
            }
        } while (changed); // Keep simplifying until no more changes are made
        return current;
    }


export NodePtr simplify(NodePtr node) {
    if (!node) return nullptr;
    return simplify_recursive(node);
}

export NodePtr parse_function_string(const std::string& input) {
    Parser parser(input);
    NodePtr root = parser.parse();
    // root->original_expression = input; // Consider if this member is still useful/maintained
    return root;
}