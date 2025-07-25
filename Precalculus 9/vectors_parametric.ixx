export module vectors_parametric;

import <vector>;
import <cmath>;
import <optional>;
import <concepts>;
import <tuple>;
export constexpr double PI = 3.14159265358979323846;

export namespace vectors {

    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    // Vector types
    export template <Arithmetic T>
        using Vec2 = std::pair<T, T>;

    export template <Arithmetic T>
        using Vec3 = std::tuple<T, T, T>;

    // --- Vector Operations ---
    export template <Arithmetic T>
        Vec2<T> add(Vec2<T> a, Vec2<T> b) {
        return { a.first + b.first, a.second + b.second };
    }

    export template <Arithmetic T>
        Vec3<T> add(Vec3<T> a, Vec3<T> b) {
        return {
            std::get<0>(a) + std::get<0>(b),
            std::get<1>(a) + std::get<1>(b),
            std::get<2>(a) + std::get<2>(b)
        };
    }

    export template <Arithmetic T>
        Vec2<T> scale(Vec2<T> v, T k) {
        return { v.first * k, v.second * k };
    }

    export template <Arithmetic T>
        Vec3<T> scale(Vec3<T> v, T k) {
        return {
            std::get<0>(v) * k,
            std::get<1>(v) * k,
            std::get<2>(v) * k
        };
    }

    export template <Arithmetic T>
        T magnitude(Vec2<T> v) {
        return std::sqrt(v.first * v.first + v.second * v.second);
    }

    export template <Arithmetic T>
        T magnitude(Vec3<T> v) {
        return std::sqrt(
            std::get<0>(v) * std::get<0>(v) +
            std::get<1>(v) * std::get<1>(v) +
            std::get<2>(v) * std::get<2>(v)
        );
    }

    // --- Dot Product and Angle ---
    export template <Arithmetic T>
        T dot(Vec2<T> a, Vec2<T> b) {
        return a.first * b.first + a.second * b.second;
    }

    export template <Arithmetic T>
        T dot(Vec3<T> a, Vec3<T> b) {
        return
            std::get<0>(a) * std::get<0>(b) +
            std::get<1>(a) * std::get<1>(b) +
            std::get<2>(a) * std::get<2>(b);
    }

    export template <Arithmetic T>
        std::optional<T> angle_between(Vec2<T> a, Vec2<T> b) {
        T magA = magnitude(a);
        T magB = magnitude(b);
        if (magA == 0 || magB == 0) return std::nullopt;

        T cos_theta = dot(a, b) / (magA * magB);
        // Clamp cos_theta to the valid range [-1, 1]
        if (cos_theta > 1.0) cos_theta = 1.0;
        else if (cos_theta < -1.0) cos_theta = -1.0;

        return std::acos(cos_theta);
    }

    // --- Vector Projection ---
    export template <Arithmetic T>
        std::optional<Vec2<T>> projection(Vec2<T> a, Vec2<T> b) {
        T mag2 = dot(b, b);
        if (mag2 == 0) return std::nullopt;
        return scale(b, dot(a, b) / mag2);
    }

    // --- Parametric Equations ---
    export template <Arithmetic T>
        Vec2<T> parametric_position(Vec2<T> initial, Vec2<T> velocity, T t) {
        return {
            initial.first + velocity.first * t,
            initial.second + velocity.second * t
        };
    }

    // --- Eliminating Parameter (only if x(t) is linear and invertible) ---
    export template <Arithmetic T>
        std::optional<std::pair<T, T>> eliminate_parameter(Vec2<T> initial, Vec2<T> velocity) {
        if (velocity.first == 0) return std::nullopt; // x constant
        T m = velocity.second / velocity.first;
        T b = initial.second - m * initial.first;
        return std::pair<T, T>{ m, b }; // y = mx + b
    }
}
