#include <iostream>
#include <vector>
#include <cmath>
#include <concepts>
#include <optional>
#include <iomanip>
#include <functional>
#include <string>
#include <numeric>
/*
3. Exponential and Logarithmic Functions
Exponential functions and their properties

Logarithms and log properties

Solving exponential and logarithmic equations

Applications of exponential/logarithmic models (growth/decay)
*/
template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Evaluate function on domain
template<Arithmetic T>
std::vector<std::optional<T>> evaluate(
    const std::function<std::optional<T>(T)>& func,
    const std::vector<T>& domain
) {
    std::vector<std::optional<T>> out;
    out.reserve(domain.size());
    for (T x : domain) {
        out.push_back(func(x));
    }
    return out;
}

// Output table
template<Arithmetic T>
void print_table(const std::vector<T>& domain, const std::vector<std::optional<T>>& range, const std::string& name) {
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\nFunction: " << name << "\n";
    std::cout << "x\tf(x)\n";
    for (size_t i = 0; i < domain.size(); ++i) {
        std::cout << domain[i] << "\t";
        if (!range[i].has_value()) std::cout << "undefined\n";
        else std::cout << range[i].value() << "\n";
    }
}

// Solve exponential equations: a^x = b → x = log_b(b)/log_b(a)
template<Arithmetic T>
std::optional<T> solve_exp_equation(T a, T b) {
    if (a <= 0 || a == 1 || b <= 0) return std::nullopt;
    return std::log(b) / std::log(a);
}

// Solve log equations: log_a(x) = b → x = a^b
template<Arithmetic T>
std::optional<T> solve_log_equation(T a, T b) {
    if (a <= 0 || a == 1) return std::nullopt;
    return std::pow(a, b);
}

// Exponential growth/decay: A(t) = A0 * e^(kt)
template<Arithmetic T>
std::function<std::optional<T>(T)> exponential_model(T A0, T k) {
    return [=](T t) -> std::optional<T> {
        return A0 * std::exp(k * t);
        };
}

int main() {
    using T = double;

    // Domain for plotting
    std::vector<T> domain;
    for (T x = -5.0; x <= 5.0; x += 0.5) domain.push_back(x);

    // Exponential functions
    auto exp_base2 = [](T x) -> std::optional<T> {
        return std::pow(2.0, x);
        };
    auto exp_base_e = [](T x) -> std::optional<T> {
        return std::exp(x);
        };
    auto decay = [](T x) -> std::optional<T> {
        return std::exp(-x);
        };

    // Logarithmic functions
    auto log_base10 = [](T x) -> std::optional<T> {
        if (x > 0) return std::log10(x); else return std::nullopt;
        };
    auto log_base_e = [](T x) -> std::optional<T> {
        if (x > 0) return std::log(x); else return std::nullopt;
        };
    auto log_base2 = [](T x) -> std::optional<T> {
        if (x > 0) return std::log(x) / std::log(2.0); else return std::nullopt;
        };

    // Evaluate and print
    print_table(domain, evaluate<T>(exp_base2, domain), "2^x");
    print_table(domain, evaluate<T>(exp_base_e, domain), "e^x");
    print_table(domain, evaluate<T>(decay, domain), "e^(-x)");

    print_table(domain, evaluate<T>(log_base10, domain), "log10(x)");
    print_table(domain, evaluate<T>(log_base_e, domain), "ln(x)");
    print_table(domain, evaluate<T>(log_base2, domain), "log2(x)");

    // Solving equations
    auto sol1 = solve_exp_equation<T>(2.0, 16.0);  // 2^x = 16 → x = 4
    auto sol2 = solve_log_equation<T>(10.0, 2.0);  // log10(x) = 2 → x = 100

    std::cout << "\nSolving equations:\n";
    std::cout << "2^x = 16 → x = " << (sol1.has_value() ? std::to_string(sol1.value()) : "undefined") << "\n";
    std::cout << "log10(x) = 2 → x = " << (sol2.has_value() ? std::to_string(sol2.value()) : "undefined") << "\n";

    // Application: exponential growth (e.g., population)
    T A0 = 1000;
    T k_growth = 0.3;  // growth rate
    auto growth_model = exponential_model<T>(A0, k_growth);
    print_table(domain, evaluate<T>(growth_model, domain), "Population Growth: A(t) = 1000 * e^(0.3t)");

    // Application: exponential decay (e.g., radioactive)
    T k_decay = -0.5;
    auto decay_model = exponential_model<T>(A0, k_decay);
    print_table(domain, evaluate<T>(decay_model, domain), "Decay: A(t) = 1000 * e^(-0.5t)");

    return 0;
}
