/*
4. Trigonometric Functions
    Angles and radians

    Unit circle definition

    Graphs of sine, cosine, tangent, etc.

    Inverse trig functions

    Trigonometric identities
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <concepts>
#include <iomanip>
#include <optional>
#include <string>
#include <functional>

constexpr double PI = 3.14159265358979323846;

template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// Degree to radian conversion
template<Arithmetic T>
constexpr T deg_to_rad(T degrees) {
    return degrees * PI / 180.0;
}

// Radian to degree conversion
template<Arithmetic T>
constexpr T rad_to_deg(T radians) {
    return radians * 180.0 / PI;
}

// Evaluate over domain
template<Arithmetic T>
std::vector<std::optional<T>> evaluate(
    const std::function<std::optional<T>(T)>& func,
    const std::vector<T>& domain
) {
    std::vector<std::optional<T>> result;
    result.reserve(domain.size());
    for (T x : domain)
        result.push_back(func(x));
    return result;
}

// Print function table
template<Arithmetic T>
void print_table(
    const std::vector<T>& domain,
    const std::vector<std::optional<T>>& range,
    const std::string& name,
    bool degrees = false
) {
    std::cout << "\n" << name << "\n";
    std::cout << (degrees ? "θ (deg)" : "x (rad)") << "\tf(x)\n";
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < domain.size(); ++i) {
        std::cout << (degrees ? rad_to_deg(domain[i]) : domain[i]) << "\t";
        if (!range[i].has_value()) std::cout << "undefined\n";
        else std::cout << range[i].value() << "\n";
    }
}

// Check identity: f(x) == g(x) within epsilon
template<Arithmetic T>
void verify_identity(
    const std::function<std::optional<T>(T)>& f,
    const std::function<std::optional<T>(T)>& g,
    const std::vector<T>& domain,
    const std::string& identity_name,
    T epsilon = 1e-6
) {
    bool holds = true;
    for (T x : domain) {
        auto fx = f(x), gx = g(x);
        if (!fx.has_value() || !gx.has_value()) continue;
        if (std::abs(fx.value() - gx.value()) > epsilon) {
            holds = false;
            break;
        }
    }
    std::cout << "Identity [" << identity_name << "] " << (holds ? "holds ✅" : "fails ❌") << "\n";
}

int main() {
    using T = double;

    // Domain: -2pi to 2pi
    std::vector<T> rad_domain;
    for (T x = -2 * PI; x <= 2 * PI; x += PI / 12.0)
        rad_domain.push_back(x);

    // Angle ↔ radian tests
    std::cout << "Angle Conversions:\n";
    std::cout << "180(degrees) = " << deg_to_rad<T>(180.0) << " rad\n";
    std::cout << "pi rad = " << rad_to_deg<T>(PI) << "(degrees)\n";

    // Trig functions
    auto sin_fn = [](T x) -> std::optional<T> { return std::sin(x); };
    auto cos_fn = [](T x) -> std::optional<T> { return std::cos(x); };
    auto tan_fn = [](T x) -> std::optional<T> {
        if (std::fmod(std::abs(x - PI / 2), PI) < 1e-6) return std::nullopt; // asymptote
        return std::tan(x);
        };

    auto arcsin_fn = [](T x) -> std::optional<T> {
        if (x < -1.0 || x > 1.0) return std::nullopt;
        return std::asin(x);
        };
    auto arccos_fn = [](T x) -> std::optional<T> {
        if (x < -1.0 || x > 1.0) return std::nullopt;
        return std::acos(x);
        };
    auto arctan_fn = [](T x) -> std::optional<T> {
        return std::atan(x);
        };

    print_table(rad_domain, evaluate<T>(sin_fn, rad_domain), "sin(x)");
    print_table(rad_domain, evaluate<T>(cos_fn, rad_domain), "cos(x)");
    print_table(rad_domain, evaluate<T>(tan_fn, rad_domain), "tan(x)");

    // Unit circle points every 30(degrees)
    std::vector<T> unit_angles;
    for (int deg = 0; deg <= 360; deg += 30)
        unit_angles.push_back(deg_to_rad<T>(deg));
    print_table(unit_angles, evaluate<T>(sin_fn, unit_angles), "sin(θ) on Unit Circle", true);
    print_table(unit_angles, evaluate<T>(cos_fn, unit_angles), "cos(θ) on Unit Circle", true);

    // Inverse trig
    std::vector<T> inv_domain;
    for (T x = -1.0; x <= 1.0; x += 0.1)
        inv_domain.push_back(x);
    print_table(inv_domain, evaluate<T>(arcsin_fn, inv_domain), "arcsin(x)");
    print_table(inv_domain, evaluate<T>(arccos_fn, inv_domain), "arccos(x)");

    std::vector<T> inv_tan_domain;
    for (T x = -5.0; x <= 5.0; x += 0.5)
        inv_tan_domain.push_back(x);
    print_table(inv_tan_domain, evaluate<T>(arctan_fn, inv_tan_domain), "arctan(x)");

    // Identity verification
    verify_identity<T>(
        [](T x) -> std::optional<T> { return std::sin(x) * std::sin(x) + std::cos(x) * std::cos(x); },
        [](T x) -> std::optional<T> { return 1.0; },
        rad_domain,
        "sin^2(x) + cos^2(x) = 1"
    );

    verify_identity<T>(
        [](T x) -> std::optional<T> { return 1.0 / std::cos(x); },
        [](T x) -> std::optional<T> {
            if (std::abs(std::cos(x)) < 1e-8) return std::nullopt;
            return 1.0f / std::cos(x);
        },
        rad_domain,
        "sec(x) = 1/cos(x)"
    );

    return 0;
}
