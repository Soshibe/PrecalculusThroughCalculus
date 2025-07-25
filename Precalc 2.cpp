/*
	2. Polynomial and Rational Functions
	Power functions

	Polynomial functions and their graphs

	Dividing polynomials

	Zeros of polynomials (including the Rational Root Theorem and Descartes’ Rule)

	Graphing rational functions

	Asymptotes and discontinuities
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <concepts>
#include <optional>
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>

template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Evaluate a polynomial at x
template<Arithmetic T>
T evaluate_polynomial(const std::vector<T>& coeffs, T x) {
    T result = 0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result = result * x + coeffs[i]; // Horner's Method
    }
    return result;
}

// Vectorized polynomial evaluation
template<Arithmetic T>
std::vector<T> evaluate_over_domain(const std::vector<T>& coeffs, const std::vector<T>& domain) {
    std::vector<T> result;
    result.reserve(domain.size());
    for (T x : domain)
        result.push_back(evaluate_polynomial(coeffs, x));
    return result;
}

// Polynomial long division: returns {quotient, remainder}
template<Arithmetic T>
std::pair<std::vector<T>, std::vector<T>> divide_polynomials(const std::vector<T>& dividend, const std::vector<T>& divisor) {
    std::vector<T> a = dividend;
    std::vector<T> b = divisor;
    std::vector<T> q;

    while (a.size() >= b.size()) {
        T coeff = a[0] / b[0];
        q.push_back(coeff);

        for (size_t i = 0; i < b.size(); ++i)
            a[i] -= coeff * b[i];

        a.erase(a.begin());
    }

    return { q, a }; // remainder is what's left of a
}

// Get rational root candidates via Rational Root Theorem
std::set<int> rational_root_candidates(const std::vector<int>& coeffs) {
    int a0 = coeffs.back();
    int an = coeffs.front();
    std::set<int> p, q;

    for (int i = 1; i <= std::abs(a0); ++i)
        if (a0 % i == 0) p.insert(i), p.insert(-i);
    for (int i = 1; i <= std::abs(an); ++i)
        if (an % i == 0) q.insert(i), q.insert(-i);

    std::set<int> roots;
    for (int num : p) {
        for (int den : q) {
            if (den == 0) continue;
            roots.insert(num / std::gcd(num, den)); // Simplify
        }
    }

    return roots;
}

// Count sign changes in polynomial coefficients (Descartes’ Rule)
int count_sign_changes(const std::vector<int>& coeffs) {
    int changes = 0;
    for (size_t i = 1; i < coeffs.size(); ++i) {
        if ((coeffs[i - 1] > 0 && coeffs[i] < 0) || (coeffs[i - 1] < 0 && coeffs[i] > 0))
            ++changes;
    }
    return changes;
}

// Rational function f(x) = P(x) / Q(x), with domain handling
template<Arithmetic T>
std::vector<std::optional<T>> evaluate_rational(
    const std::vector<T>& numerator,
    const std::vector<T>& denominator,
    const std::vector<T>& domain
) {
    std::vector<std::optional<T>> result;
    result.reserve(domain.size());

    for (T x : domain) {
        T num = evaluate_polynomial(numerator, x);
        T den = evaluate_polynomial(denominator, x);
        if (std::abs(den) < 1e-8)
            result.push_back(std::nullopt);
        else
            result.push_back(num / den);
    }
    return result;
}

// Pretty print
template<Arithmetic T>
void print_function_table(
    const std::vector<T>& domain,
    const std::vector<std::optional<T>>& range,
    const std::string& name
) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n" << name << "\n";
    std::cout << "x\tf(x)\n";
    for (size_t i = 0; i < domain.size(); ++i) {
        std::cout << domain[i] << "\t";
        if (!range[i].has_value())
            std::cout << "undefined\n";
        else
            std::cout << range[i].value() << "\n";
    }
}

int main() {
    using T = double;

    // Domain for plotting
    std::vector<T> domain;
    for (T x = -10; x <= 10; x += 1.0)
        domain.push_back(x);

    // Power function: f(x) = x^5
    std::vector<T> power_coeffs = { 1, 0, 0, 0, 0, 0 }; // x^5
    auto power_range = evaluate_over_domain<T>(power_coeffs, domain);
    std::vector<std::optional<T>> power_opt(power_range.begin(), power_range.end());
    print_function_table(domain, power_opt, "Power Function: x^5");

    // Polynomial: f(x) = 2x^3 - 3x^2 - 11x + 6
    std::vector<T> poly = { 2, -3, -11, 6 };
    auto poly_range = evaluate_over_domain<T>(poly, domain);
    std::vector<std::optional<T>> poly_opt(poly_range.begin(), poly_range.end());
    print_function_table(domain, poly_opt, "Polynomial f(x) = 2x^3 - 3x^2 - 11x + 6");

    // Polynomial division: divide f(x) by (x - 1)
    std::vector<T> divisor = { 1, -1 }; // x - 1
    auto [quotient, remainder] = divide_polynomials<T>(poly, divisor);
    std::cout << "\nQuotient: ";
    for (auto c : quotient) std::cout << c << " ";
    std::cout << "\nRemainder: ";
    for (auto c : remainder) std::cout << c << " ";
    std::cout << "\n";

    // Rational Root Theorem (integer coeffs only)
    std::vector<int> poly_int = { 2, -3, -11, 6 };
    auto candidates = rational_root_candidates(poly_int);
    std::cout << "\nRational Root Candidates:\n";
    for (auto r : candidates) std::cout << r << " ";
    std::cout << "\n";

    // Descartes’ Rule of Signs
    int sign_changes = count_sign_changes(poly_int);
    std::cout << "Descartes' Sign Changes: " << sign_changes << "\n";

    // Rational function: f(x) = (x^2 - 1) / (x - 1)
    std::vector<T> num = { 1, 0, -1 };  // x^2 - 1
    std::vector<T> den = { 1, -1 };     // x - 1
    auto rational = evaluate_rational<T>(num, den, domain);
    print_function_table(domain, rational, "Rational Function: (x^2 - 1)/(x - 1)");

    return 0;
}
