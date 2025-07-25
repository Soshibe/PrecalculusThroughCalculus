export module trig_applications;

import <cmath>;
import <vector>;
import <optional>;
import <concepts>;
import <functional>;
import <string>;
import <iostream>;
import <iomanip>;

#if __cplusplus >= 202002L
#include <numbers>
export constexpr double PI = std::numbers::pi_v<double>;
#else
export constexpr double PI = 3.14159265358979323846;
#endif

export constexpr double EPS = 1e-6;

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

// --- Law of Sines ---

// Computes missing side c given a, A, and C (a / sin A = c / sin C)
export template <Arithmetic T>
std::optional<T> law_of_sines_side(T a, T A_rad, T C_rad) {
    if (a <= 0 || A_rad <= 0 || C_rad <= 0 || A_rad >= PI || C_rad >= PI)
        return std::nullopt;
    T sinA = std::sin(A_rad);
    T sinC = std::sin(C_rad);
    if (sinA < EPS) return std::nullopt;
    return a * sinC / sinA;
}

// Computes missing angle C given sides a, c and angle A using Law of Sines
export template <Arithmetic T>
std::optional<T> law_of_sines_angle(T a, T c, T A_rad) {
    if (a <= 0 || c <= 0 || A_rad <= 0 || A_rad >= PI)
        return std::nullopt;
    T sinA = std::sin(A_rad);
    T ratio = c * sinA / a;
    if (ratio < -1.0 - EPS || ratio > 1.0 + EPS)
        return std::nullopt;
    if (ratio > 1.0) ratio = 1.0;
    else if (ratio < -1.0) ratio = -1.0;
    return std::asin(ratio);
}

// --- Law of Cosines ---

// Computes side c from sides a, b and angle C (in radians)
export template <Arithmetic T>
std::optional<T> law_of_cosines_side(T a, T b, T C_rad) {
    if (a <= 0 || b <= 0 || C_rad < 0 || C_rad >= PI)
        return std::nullopt;
    T c_sq = a * a + b * b - 2 * a * b * std::cos(C_rad);
    if (c_sq < 0) return std::nullopt;
    return std::sqrt(c_sq);
}

// Computes angle C from sides a, b, c using Law of Cosines
export template <Arithmetic T>
std::optional<T> law_of_cosines_angle(T a, T b, T c) {
    if (a <= 0 || b <= 0 || c <= 0)
        return std::nullopt;
    T numerator = a * a + b * b - c * c;
    T denominator = 2 * a * b;
    if (denominator < EPS) return std::nullopt;
    T val = numerator / denominator;
    if (val < -1.0 - EPS || val > 1.0 + EPS) return std::nullopt;
    if (val > 1.0) val = 1.0;
    else if (val < -1.0) val = -1.0;
    return std::acos(val);
}

// --- Area of Triangle ---

// Using 1/2 * a * b * sin C
export template <Arithmetic T>
std::optional<T> triangle_area_sas(T a, T b, T C_rad) {
    if (a <= 0 || b <= 0 || C_rad <= 0 || C_rad >= PI)
        return std::nullopt;
    return T(0.5) * a * b * std::sin(C_rad);
}

// Using Heron's formula given sides a,b,c
export template <Arithmetic T>
std::optional<T> triangle_area_heron(T a, T b, T c) {
    if (a <= 0 || b <= 0 || c <= 0) return std::nullopt;
    T s = (a + b + c) / 2.0;
    T val = s * (s - a) * (s - b) * (s - c);
    if (val < 0) return std::nullopt;
    return std::sqrt(val);
}

// --- Triangle Solver ---

// Structure to hold solved triangle data
export template <Arithmetic T>
struct Triangle {
    T a, b, c;   // sides
    T A, B, C;   // angles in radians
};

// Solve triangle given SSS (all sides)
export template <Arithmetic T>
std::optional<Triangle<T>> solve_triangle_sss(T a, T b, T c) {
    if (a <= 0 || b <= 0 || c <= 0) return std::nullopt;
    auto A = law_of_cosines_angle(a, b, c);
    auto B = law_of_cosines_angle(b, c, a);
    auto C = law_of_cosines_angle(c, a, b);
    if (!A.has_value() || !B.has_value() || !C.has_value()) return std::nullopt;
    return Triangle<T>{a, b, c, A.value(), B.value(), C.value()};
}

// Solve triangle given SAS (two sides and included angle)
export template <Arithmetic T>
std::optional<Triangle<T>> solve_triangle_sas(T a, T C_rad, T b) {
    auto c_opt = law_of_cosines_side(a, b, C_rad);
    if (!c_opt.has_value()) return std::nullopt;
    auto A_opt = law_of_cosines_angle(b, c_opt.value(), a);
    auto B_opt = law_of_cosines_angle(a, c_opt.value(), b);
    if (!A_opt.has_value() || !B_opt.has_value()) return std::nullopt;
    return Triangle<T>{a, b, c_opt.value(), A_opt.value(), B_opt.value(), C_rad};
}

// Solve triangle given ASA or AAS
// input: two angles and one side opposite one angle
export template <Arithmetic T>
std::optional<Triangle<T>> solve_triangle_asa_aas(T A_rad, T B_rad, T a) {
    if (A_rad <= 0 || B_rad <= 0 || a <= 0) return std::nullopt;
    T C_rad = PI - A_rad - B_rad;
    if (C_rad <= 0) return std::nullopt;
    T b = law_of_sines_side(a, A_rad, B_rad).value_or(-1);
    T c = law_of_sines_side(a, A_rad, C_rad).value_or(-1);
    if (b <= 0 || c <= 0) return std::nullopt;
    return Triangle<T>{a, b, c, A_rad, B_rad, C_rad};
}

// Solve triangle SSA (ambiguous case) - returns zero, one, or two solutions
export template <Arithmetic T>
std::vector<Triangle<T>> solve_triangle_ssa(T a, T b, T A_rad) {
    std::vector<Triangle<T>> solutions;
    if (a <= 0 || b <= 0 || A_rad <= 0 || A_rad >= PI) return solutions;

    T sinA = std::sin(A_rad);
    T ratio = b * sinA / a;
    if (ratio < -1.0 - EPS || ratio > 1.0 + EPS) return solutions; // no solution
    if (ratio > 1.0) ratio = 1.0;
    else if (ratio < -1.0) ratio = -1.0;

    T B1 = std::asin(ratio);
    T B2 = PI - B1;

    T C1 = PI - A_rad - B1;
    if (C1 > 0) {
        T c1 = law_of_sines_side(a, A_rad, C1).value_or(-1);
        if (c1 > 0)
            solutions.push_back({ a, b, c1, A_rad, B1, C1 });
    }
    if (std::abs(B2 - B1) > EPS) { // check distinct solution
        T C2 = PI - A_rad - B2;
        if (C2 > 0) {
            T c2 = law_of_sines_side(a, A_rad, C2).value_or(-1);
            if (c2 > 0)
                solutions.push_back({ a, b, c2, A_rad, B2, C2 });
        }
    }
    return solutions;
}

// --- Trigonometric Form of Complex Numbers ---

export template <Arithmetic T>
struct ComplexPolar {
    T r;    // modulus
    T theta; // argument in radians
};

export template <Arithmetic T>
ComplexPolar<T> to_polar(T x, T y) {
    return ComplexPolar<T>{
        static_cast<T>(std::hypot(x, y)),
            static_cast<T>(std::atan2(y, x))
    };
}

export template <Arithmetic T>
std::pair<T, T> to_cartesian(ComplexPolar<T> p) {
    return { p.r * std::cos(p.theta), p.r * std::sin(p.theta) };
}

export template <Arithmetic A, Arithmetic B>
ComplexPolar<std::common_type_t<A, B>> multiply_polar(const ComplexPolar<A>& a, const ComplexPolar<B>& b) {
    using Result = std::common_type_t<A, B>;
    Result r = static_cast<Result>(a.r) * static_cast<Result>(b.r);
    Result theta = std::fmod(static_cast<Result>(a.theta) + static_cast<Result>(b.theta), 2 * PI);
    if (theta < 0) theta += 2 * PI;
    return ComplexPolar<Result>{r, theta};
}


export template <Arithmetic T>
ComplexPolar<T> divide_polar(const ComplexPolar<T>& a, const ComplexPolar<T>& b) {
    if (b.r == 0) return ComplexPolar<T>{0, 0}; // or throw
    T angle = std::fmod(a.theta - b.theta, 2 * PI);
    if (angle < 0) angle += 2 * PI;
    return ComplexPolar<T>{a.r / b.r, angle};
}

// --- Polar Coordinates ---

// Converts polar coordinates to Cartesian
export template <Arithmetic T>
std::pair<T, T> polar_to_cartesian(T r, T theta) {
    return { r * std::cos(theta), r * std::sin(theta) };
}

// Converts Cartesian coordinates to polar
export template <Arithmetic T>
std::pair<T, T> cartesian_to_polar(T x, T y) {
    return { std::hypot(x, y), std::atan2(y, x) };
}

