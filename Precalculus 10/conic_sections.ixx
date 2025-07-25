export module conic_sections;

import <cmath>;
import <optional>;
import <tuple>;
import <string>; // Not strictly used in current functions, but kept if intended for string conversions
import <concepts>;
// import <numbers>; // Uncomment and use for std::numbers::pi if available and desired

export constexpr double PI = 3.14159265358979323846; // Or use std::numbers::pi

export namespace conics {

    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    // --- Discriminant classification ---
    export enum class ConicType { Circle, Parabola, Ellipse, Hyperbola, Degenerate, Unknown };

    // Classifies the conic section based on its quadratic terms (shape),
    // from the general equation Ax² + Bxy + Cy² + Dx + Ey + F = 0.
    // Note: This classification determines the *type* of conic (e.g., parabolic shape),
    // but does not distinguish between degenerate and non-degenerate forms
    // (e.g., a point ellipse vs. a true ellipse). Full degeneracy checks require
    // analyzing all coefficients (A-F).
    export template <Arithmetic T>
        ConicType classify_general_conic(T A, T B, T C) {
        T D = B * B - 4 * A * C;

        // Handle cases where all quadratic coefficients are zero (linear equation or constant)
        if (std::abs(A) < 1e-9 && std::abs(B) < 1e-9 && std::abs(C) < 1e-9) return ConicType::Degenerate; // Represents a line or a constant equation

        if (D < 0) {
            // Elliptic type
            // Use a small epsilon for floating point comparison if T is float/double
            if (std::abs(A - C) < 1e-9 && std::abs(B) < 1e-9) return ConicType::Circle;
            return ConicType::Ellipse;
        }
        if (D == 0) return ConicType::Parabola;

        // D > 0
        return ConicType::Hyperbola;
    }

    // --- Parabolas ---
    // Returns the vertex (h, k) of a parabola given in standard form y = ax^2 + bx + c
    export template <Arithmetic T>
        std::pair<T, T> vertex_from_standard(T a, T b, T c) {
        T h = -b / (2 * a);
        T k = a * h * h + b * h + c;
        return { h, k };
    }

    // Returns the coefficients (A, B, C) for the standard quadratic form
    // y = Ax^2 + Bx + C, given its vertex (h, k) and leading coefficient 'a'.
    export template <Arithmetic T>
        std::tuple<T, T, T> parabola_from_vertex(T h, T k, T a) {
        // Corresponds to y = a(x - h)^2 + k
        // y = a(x^2 - 2hx + h^2) + k
        // y = ax^2 - 2ahx + ah^2 + k
        return std::tuple{ a, -2 * a * h, a * h * h + k };
    }

    // --- Ellipse ---
    export template <Arithmetic T>
        struct Ellipse {
        T h, k; // Center (h, k)
        T a, b; // Semi-major and semi-minor axis lengths (or vice-versa depending on orientation)
        bool horizontal_major_axis; // true if major axis is along x, false if along y

        // Constructor for convenience
        constexpr Ellipse(T h_val = 0, T k_val = 0, T a_val = 0, T b_val = 0, bool horz = true)
            : h(h_val), k(k_val), a(a_val), b(b_val), horizontal_major_axis(horz) {
        }
    };

    // Helper to create an Ellipse struct (optional if constructors are used directly)
    export template <Arithmetic T>
        Ellipse<T> create_ellipse(T h, T k, T a, T b, bool horizontal) {
        return Ellipse<T>{ h, k, a, b, horizontal };
    }

    // --- Hyperbola ---
    export template <Arithmetic T>
        struct Hyperbola {
        T h, k; // Center (h, k)
        T a, b; // a is length from center to vertex, b is length from center to co-vertex (defines asymptotes)
        bool horizontal_transverse_axis; // true if transverse axis is along x, false if along y

        // Constructor for convenience
        constexpr Hyperbola(T h_val = 0, T k_val = 0, T a_val = 0, T b_val = 0, bool horz = true)
            : h(h_val), k(k_val), a(a_val), b(b_val), horizontal_transverse_axis(horz) {
        }
    };

    // Helper to create a Hyperbola struct (optional if constructors are used directly)
    export template <Arithmetic T>
        Hyperbola<T> create_hyperbola(T h, T k, T a, T b, bool horizontal) {
        return Hyperbola<T>{ h, k, a, b, horizontal };
    }

    // --- Rotated Conics (based on Ax² + Bxy + Cy² + Dx + Ey + F = 0) ---
    // Returns the angle of rotation (theta in radians) needed to eliminate the Bxy term.
    // tan(2*theta) = B / (A - C)
    export template <Arithmetic T>
        T rotation_angle(T A, T B, T C) {
        if (std::abs(B) < 1e-9) return 0; // No Bxy term, no rotation needed

        // std::atan2 handles (A-C) == 0 correctly
        return 0.5 * std::atan2(B, A - C); // in radians
    }

    // --- General form conversion helper (not a function, just a struct) ---
    export template <Arithmetic T>
        struct GeneralConic {
        T A, B, C, D, E, F;
    };

    // Convenience function to classify using the GeneralConic struct
    export template <Arithmetic T>
        ConicType identify_conic(const GeneralConic<T>& g) {
        return classify_general_conic(g.A, g.B, g.C);
    }

    // --- Focus/Directrix (for parabola y = ax^2 + bx + c) ---
    // Returns (focus_x, focus_y, directrix_y_or_x) for a parabola.
    // This is for y = ax^2 + bx + c form (vertical parabola).
    export template <Arithmetic T>
        std::tuple<T, T, T> parabola_focus_directrix(T a, T b, T c) {
        T h = -b / (2 * a);
        T k = c - (b * b) / (4 * a); // y-coordinate of vertex

        // For y = a(x-h)^2 + k, the focal length is p = 1/(4a)
        T p = 1 / (4 * a);

        // Focus is (h, k + p)
        // Directrix is y = k - p
        return { h, k + p, k - p }; // (focus x, focus y, directrix value)
    }

    // --- Eccentricity ---
    // Calculates eccentricity of an ellipse given its struct.
    // e = sqrt(1 - (minor_axis_sq / major_axis_sq))
    export template <Arithmetic T>
        std::optional<T> eccentricity_ellipse(const Ellipse<T>& el) {
        T major_sq, minor_sq;
        if (el.horizontal_major_axis) { // 'a' is semi-major, 'b' is semi-minor
            major_sq = el.a * el.a;
            minor_sq = el.b * el.b;
        }
        else { // 'b' is semi-major, 'a' is semi-minor
            major_sq = el.b * el.b;
            minor_sq = el.a * el.a;
        }

        // Handle potential division by zero or invalid input (e.g., a=0 or b=0 creating a point/line)
        if (major_sq <= 0 || minor_sq > major_sq) return std::nullopt; // Ensure major_sq is positive and minor_sq is not greater

        return std::sqrt(1 - minor_sq / major_sq);
    }

    // Calculates eccentricity of a hyperbola given its struct.
    // e = sqrt(1 + (b^2 / a^2)), where 'a' is semi-transverse and 'b' is semi-conjugate.
    export template <Arithmetic T>
        std::optional<T> eccentricity_hyperbola(const Hyperbola<T>& hyp) {
        T a_sq = hyp.a * hyp.a; // Always semi-transverse axis squared
        T b_sq = hyp.b * hyp.b; // Always semi-conjugate axis squared

        // Handle potential division by zero or invalid input (e.g., a=0)
        if (a_sq <= 0) return std::nullopt;

        return std::sqrt(1 + b_sq / a_sq);
    }
}