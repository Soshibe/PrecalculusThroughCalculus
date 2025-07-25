#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <concepts>
#include <iomanip>
#include <string>
#include <limits>
#include <optional>
/*
    1. Functions and Their Graphs
    Basics of functions

    Domain and range

    Function notation

    Types of functions (linear, quadratic, etc.)

    Transformations of functions

    Combining functions (addition, subtraction, composition)

    Inverse functions
*/

template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Evaluate a function over a domain
template<Arithmetic T>
std::vector<T> evaluate_function(const std::function<std::optional<T>(T)>& func, const std::vector<T>& domain) {
    std::vector<T> result;
    result.reserve(domain.size());
    for (T x : domain) {
        auto fx = func(x);
        result.push_back(fx.value_or(std::numeric_limits<T>::quiet_NaN()));
    }
    return result;
}

// Pretty print
template<Arithmetic T>
void print_function_table(const std::vector<T>& domain, const std::vector<T>& range, const std::string& name) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nFunction: " << name << "\n";
    std::cout << "x\tf(x)\n";
    for (size_t i = 0; i < domain.size(); ++i) {
        std::cout << domain[i] << "\t";
        if (std::isnan(range[i])) std::cout << "undefined\n";
        else std::cout << range[i] << "\n";
    }
}

// Transformations: shift, scale, reflect
template<Arithmetic T>
std::function<std::optional<T>(T)> transform_function(
    const std::function<std::optional<T>(T)>& f,
    T vertical_shift = 0,
    T horizontal_shift = 0,
    T vertical_scale = 1,
    bool reflect_x = false,
    bool reflect_y = false
) {
    return [=](T x) -> std::optional<T> {
        T x_mod = reflect_y ? -x : x;
        auto y = f(x_mod - horizontal_shift);
        if (!y.has_value()) return std::nullopt;
        T y_val = vertical_scale * y.value();
        if (reflect_x) y_val *= -1;
        return y_val + vertical_shift;
        };
}

// Function composition f(g(x))
template<Arithmetic T>
std::function<std::optional<T>(T)> compose(
    const std::function<std::optional<T>(T)>& f,
    const std::function<std::optional<T>(T)>& g
) {
    return [=](T x) -> std::optional<T> {
        auto gx = g(x);
        if (!gx.has_value()) return std::nullopt;
        return f(gx.value());
        };
}

// Check if a function is even or odd
template<Arithmetic T>
void check_symmetry(const std::function<std::optional<T>(T)>& f, const std::vector<T>& domain) {
    bool even = true, odd = true;
    for (T x : domain) {
        auto f_x = f(x);
        auto f_negx = f(-x);
        if (!f_x.has_value() || !f_negx.has_value()) continue;
        if (std::abs(f_x.value() - f_negx.value()) > 1e-5) even = false;
        if (std::abs(f_x.value() + f_negx.value()) > 1e-5) odd = false;
    }
    std::cout << "Symmetry: ";
    if (even) std::cout << "Even\n";
    else if (odd) std::cout << "Odd\n";
    else std::cout << "Neither\n";
}

// Inverse check: f⁻¹(f(x)) ~= x
template<Arithmetic T>
void check_inverse(
    const std::function<std::optional<T>(T)>& f,
    const std::function<std::optional<T>(T)>& f_inv,
    const std::vector<T>& domain
) {
    bool valid = true;
    for (T x : domain) {
        auto fx = f(x);
        if (!fx.has_value()) continue;
        auto inv = f_inv(fx.value());
        if (!inv.has_value() || std::abs(inv.value() - x) > 1e-5) {
            valid = false;
            break;
        }
    }
    std::cout << "Inverse Function: " << (valid ? "Valid\n" : "Invalid\n");
}

int main() {
    using T = double;
    std::vector<T> domain;
    for (T x = -10; x <= 10; x += 1.0) domain.push_back(x);

    // Basic function
    auto linear = [](T x) -> std::optional<T> { return 2 * x + 1; };

    // Piecewise function
    auto piecewise = [](T x) -> std::optional<T> {
        if (x < 0) return x * x;
        else if (x == 0) return 0;
        else return x + 1;
        };

    // Function with undefined points (e.g., reciprocal)
    auto reciprocal = [](T x) -> std::optional<T> {
        if (x == 0) return std::nullopt;
        return 1 / x;
        };

    // Inverse of linear
    auto inverse_linear = [](T y) -> std::optional<T> {
        return (y - 1) / 2;
        };

    // Transform linear function
    auto transformed = transform_function<T>(linear, 3.0, -2.0, 0.5, true, false);

    // Compose: f(g(x)) = linear(piecewise(x))
    auto composed = compose<T>(linear, piecewise);

    // Evaluate all
    auto eval_and_print = [&](const std::function<std::optional<T>(T)>& f, const std::string& name) {
        auto range = evaluate_function<T>(f, domain);
        print_function_table<T>(domain, range, name);
        check_symmetry<T>(f, domain);
        };

    eval_and_print(linear, "linear");
    eval_and_print(piecewise, "piecewise");
    eval_and_print(reciprocal, "reciprocal");
    eval_and_print(transformed, "transformed linear");
    eval_and_print(composed, "composed f(g(x))");

    check_inverse<T>(linear, inverse_linear, domain);

    return 0;
}
