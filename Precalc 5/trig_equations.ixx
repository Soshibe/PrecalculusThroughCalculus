export module trig_equations;

// Standard library imports for the module
import <cmath>;        // For trigonometric functions (sin, cos, tan, asin, acos, atan, fmod, sqrt, pow)
import <vector>;       // For std::vector
import <optional>;     // For std::optional to handle undefined results
import <concepts>;     // For C++20 concepts (Arithmetic)
import <functional>;   // For std::function (to pass functions as arguments)
import <string>;       // For std::string
import <iostream>;     // For std::cout, std::fixed, std::setprecision
import <iomanip>;      // For std::fixed, std::setprecision

// For std::numbers::pi_v in C++20 and later
#if __cplusplus >= 202002L
#include <numbers>
export constexpr double PI = std::numbers::pi_v<double>;
#else
export constexpr double PI = 3.14159265358979323846; // Fallback for older C++ standards
#endif

export constexpr double EPS = 1e-6; // Epsilon for floating-point comparisons

// C++20 concept to restrict template types to arithmetic types
export template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// --- Helper Functions for Angle Conversions ---

// Converts degrees to radians
export template <Arithmetic T>
T deg_to_rad(T deg) {
    return deg * PI / 180.0;
}

// Converts radians to degrees
export template <Arithmetic T>
T rad_to_deg(T rad) {
    return rad * 180.0 / PI;
}

// --- Function Evaluation ---

// Evaluates a given function `f` over a `domain` and returns results as std::optional
export template <Arithmetic T>
std::vector<std::optional<T>> evaluate(const std::function<std::optional<T>(T)>& f, const std::vector<T>& domain) {
    std::vector<std::optional<T>> result;
    result.reserve(domain.size());
    for (T x : domain)
        result.push_back(f(x));
    return result;
}

// --- Robust Trigonometric Functions (handling undefined points) ---

// Robust tangent function that returns std::nullopt for asymptotes (odd multiples of PI/2)
export template <Arithmetic T>
std::optional<T> tan_robust(T x, T epsilon = EPS) {
    // std::fmod(x, PI) gives result in (-PI, PI).
    // We care about values near PI/2 and -PI/2.
    T val_mod_pi = std::fmod(x, PI);

    // Check if near asymptotes: PI/2 + n*PI
    // Check points near PI/2 (positive direction) or -PI/2 (negative direction)
    // Note: fmod result can be negative, so abs(val_mod_pi) might be close to PI/2.
    if (std::abs(val_mod_pi - PI / 2.0) < epsilon ||
        std::abs(val_mod_pi + PI / 2.0) < epsilon) {
        return std::nullopt; // Undefined (vertical asymptote)
    }
    return std::tan(x);
}

// Robust inverse trigonometric functions handling domain restrictions
export template <Arithmetic T>
std::optional<T> arcsin(T x, T epsilon = EPS) {
    if (x < -1.0 - epsilon || x > 1.0 + epsilon)
        return std::nullopt; // Outside domain [-1, 1]
    return std::asin(x);
}

export template <Arithmetic T>
std::optional<T> arccos(T x, T epsilon = EPS) {
    if (x < -1.0 - epsilon || x > 1.0 + epsilon)
        return std::nullopt; // Outside domain [-1, 1]
    return std::acos(x);
}

export template <Arithmetic T>
std::optional<T> arctan(T x) {
    return std::atan(x); // Defined for all real x
}

export template <Arithmetic T>
std::optional<T> arccot(T x, T epsilon = EPS) {
    if (std::abs(x) < epsilon) { // Special case for arccot(0)
        return PI / 2.0; // arccot(0) = PI/2
    }
    return std::atan(1.0 / x);
}

export template <Arithmetic T>
std::optional<T> arcsec(T x, T epsilon = EPS) {
    if (std::abs(x) < 1.0 - epsilon) // Outside domain (-inf, -1] U [1, inf)
        return std::nullopt;
    return std::acos(1.0 / x);
}

export template <Arithmetic T>
std::optional<T> arccsc(T x, T epsilon = EPS) {
    if (std::abs(x) < 1.0 - epsilon) // Outside domain (-inf, -1] U [1, inf)
        return std::nullopt;
    return std::asin(1.0 / x);
}

// --- Proving Trig Identities ---

// Structure to define a trigonometric identity
export template <Arithmetic T>
struct TrigIdentity {
    std::string name;                                // Name/description of the identity
    std::function<std::optional<T>(T)> lhs;         // Left-hand side of the identity
    std::function<std::optional<T>(T)> rhs;         // Right-hand side of the identity
};

// Verifies a single trigonometric identity numerically over a given domain
export template <Arithmetic T>
void verify_identity(const std::function<std::optional<T>(T)>& f,
    const std::function<std::optional<T>(T)>& g,
    const std::vector<T>& domain,
    const std::string& name,
    T epsilon = EPS) {
    bool holds = true;
    for (T x : domain) {
        auto a = f(x), b = g(x);
        // If one side is defined and the other is not, the identity fails.
        if (a.has_value() != b.has_value()) {
            holds = false;
            break;
        }
        // If both are defined but differ significantly, the identity fails.
        if (a.has_value() && std::abs(a.value() - b.value()) > epsilon) {
            holds = false;
            break;
        }
        // If both are undefined, or both are defined and match, it holds for this point.
    }
    std::cout << "Identity [" << name << "] " << (holds ? "holds" : "fails") << "\n";
}

// Verifies a vector of trigonometric identities
export template <Arithmetic T>
void verify_all_identities(const std::vector<TrigIdentity<T>>& identities,
    const std::vector<T>& domain,
    T epsilon = EPS) {
    for (const auto& identity : identities) {
        verify_identity<T>(identity.lhs, identity.rhs, domain, identity.name, epsilon);
    }
}

// --- Solving Trigonometric Equations ---

// Provides the general symbolic solution for basic trigonometric equations (sin(x)=k, cos(x)=k, tan(x)=k)
export template <Arithmetic T>
void solve_trig_equation_symbolic(const std::string& func_name, T target) {
    std::cout << "\n[Symbolic General Solution for " << func_name << "(x) = " << target << "]\n";
    std::cout << std::fixed << std::setprecision(6); // Set precision for output

    if (func_name == "sin") {
        if (std::abs(target) > 1.0 + EPS) { // Check if target is outside [-1, 1] range (with epsilon tolerance)
            std::cout << "No real solution: |target| > 1\n";
            return;
        }
        T alpha = std::asin(target); // Principal value in [-PI/2, PI/2]
        std::cout << "x = n*pi + (-1)^n*" << alpha << " rad  ~= " << rad_to_deg(alpha) << "(degrees)\n";
    }
    else if (func_name == "cos") {
        if (std::abs(target) > 1.0 + EPS) { // Check if target is outside [-1, 1] range
            std::cout << "No real solution: |target| > 1\n";
            return;
        }
        T alpha = std::acos(target); // Principal value in [0, PI]
        std::cout << "x = +/-" << alpha << " + 2*n*pi  ~= +/-" << rad_to_deg(alpha) << "(degrees) + 360(degrees)*n\n";
    }
    else if (func_name == "tan") {
        T alpha = std::atan(target); // Principal value in (-PI/2, PI/2)
        std::cout << "x = " << alpha << " + n*pi  ~= " << rad_to_deg(alpha) << "(degrees) + 180(degrees)*n\n";
    }
    else {
        std::cout << "Unsupported function for symbolic solution.\n";
    }
}

// Numerically finds solutions for a trigonometric equation `f(x) = target` within a given range
export template <typename T>
std::vector<T> solve_trig_equation_numeric(const std::function<std::optional<T>(T)>& f,
    T target,
    T xmin, T xmax,
    T step = PI / 1800.0f, // Smaller step for more precision, ~1.875 degrees
    T epsilon = EPS) {
    std::vector<T> roots;
    for (T x = xmin; x <= xmax; x += step) {
        auto y = f(x);
        if (!y.has_value()) continue; // Skip if function is undefined at this point

        if (std::abs(y.value() - target) < epsilon) { // Check if function value is close to target
            // Deduplicate solutions: check if the current root is significantly different
            // from previously found roots to avoid capturing duplicates due to step size.
            bool unique = true;
            for (auto r : roots) {
                if (std::abs(r - x) < step / 2.0) { // If within half a step of an existing root
                    unique = false;
                    break;
                }
            }
            if (unique) roots.push_back(x);
        }
    }
    return roots;
}

// Helper to pretty-print numeric solutions (in radians and degrees)
export template <typename T>
void print_numeric_solutions(const std::vector<T>& roots) {
    std::cout << std::fixed << std::setprecision(6);
    for (T x : roots) {
        std::cout << "x ~= " << x << " rad  ~=  " << rad_to_deg(x) << "(degrees)\n";
    }
}

// The 'solve_trig_equation' function that takes a pre-generated domain.
// Keeping it for completeness if there's a specific use case, otherwise
// 'solve_trig_equation_numeric' is generally more versatile.
export template <Arithmetic T>
void solve_trig_equation(const std::function<std::optional<T>(T)>& f,
    const std::vector<T>& domain,
    const std::string& name,
    T target,
    T epsilon = EPS) {
    std::cout << "\nSolutions for " << name << " = " << target << "\n";
    for (T x : domain) {
        auto val = f(x);
        if (val.has_value() && std::abs(val.value() - target) < epsilon)
            std::cout << "x ~= " << x << " rad ~= " << rad_to_deg(x) << "(degrees)\n";
    }
}