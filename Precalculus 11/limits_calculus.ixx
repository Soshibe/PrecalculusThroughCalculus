export module limits_calculus;

#include <cmath>      // For std::abs, std::isnan, std::isinf, std::sqrt (though not used directly here)
#include <optional>   // For std::optional
#include <concepts>   // For std::is_arithmetic_v
#include <limits>     // For std::numeric_limits (not directly used in the final version, but useful for understanding epsilon)

export namespace calculus {

    /**
     * @brief Concept for arithmetic types (integral or floating-point).
     * @tparam T The type to check.
     */
    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    /**
     * @brief Structure to hold the results of a numerical limit approximation.
     * @tparam T The arithmetic type of the function's output.
     */
    export template <Arithmetic T>
        struct LimitResult {
        std::optional<T> left_limit;  ///< The approximated limit from the left side. std::nullopt if not computable.
        std::optional<T> right_limit; ///< The approximated limit from the right side. std::nullopt if not computable.
        bool limits_agree;            ///< True if both limits are present and agree within the specified tolerance.

        /**
         * @brief Checks if a single, two-sided limit exists and is finite.
         * This is true if both left and right limits are present and agree.
         * @return True if a finite, two-sided limit is approximated to exist, false otherwise.
         */
        bool has_two_sided_limit() const {
            return left_limit.has_value() && right_limit.has_value() && limits_agree;
        }

        /**
         * @brief Returns the value of the two-sided limit if it exists and agrees.
         * Returns the average of left and right limits for a better approximation.
         * @return An optional containing the two-sided limit, or std::nullopt if it doesn't exist or agree.
         */
        std::optional<T> get_two_sided_limit() const {
            if (has_two_sided_limit()) {
                return (left_limit.value() + right_limit.value()) / static_cast<T>(2.0);
            }
            return std::nullopt;
        }
    };

    /**
     * @brief Numerically approximates the limit of a function f(x) as x approaches x0
     * from the left and right sides.
     * @tparam T The arithmetic type of the function's input and output.
     * @tparam Func The type of the callable function (e.g., lambda, std::function).
     * @param f The function for which to approximate the limit.
     * @param x0 The point at which to approximate the limit.
     * @param delta A small positive value representing the step size away from x0.
     * A default of 1e-6 is chosen for common floating-point precision.
     * @param tolerance The maximum allowed absolute difference between the left and right
     * limit approximations for them to be considered "agreeing".
     * A default of 1e-5 is chosen (10x delta, a common heuristic).
     * @return A LimitResult structure containing the approximated left and right limits,
     * and a boolean indicating if they agree within the tolerance.
     * Returns std::nullopt for a limit value if the function evaluates to NaN or Inf
     * at the probed points.
     */
    export template <Arithmetic T, typename Func>
        LimitResult<T> approximate_limit(Func f, T x0, T delta = static_cast<T>(1e-6), T tolerance = static_cast<T>(1e-5)) {
        if (delta <= 0) {
            // Delta must be positive to probe points around x0.
            return { std::nullopt, std::nullopt, false };
        }
        if (tolerance < 0) { // Tolerance should be non-negative
            tolerance = std::abs(tolerance);
        }

        T left_x = x0 - delta;
        T right_x = x0 + delta;

        std::optional<T> left_val_opt;
        std::optional<T> right_val_opt;

        // Attempt to evaluate f(left_x)
        T temp_left_val = f(left_x);
        if (!std::isnan(temp_left_val) && !std::isinf(temp_left_val)) {
            left_val_opt = temp_left_val;
        }

        // Attempt to evaluate f(right_x)
        T temp_right_val = f(right_x);
        if (!std::isnan(temp_right_val) && !std::isinf(temp_right_val)) {
            right_val_opt = temp_right_val;
        }

        bool agree = false;
        if (left_val_opt.has_value() && right_val_opt.has_value()) {
            agree = std::abs(left_val_opt.value() - right_val_opt.value()) < tolerance;
        }

        return { left_val_opt, right_val_opt, agree };
    }

    /**
     * @brief Numerically approximates the derivative of a function f(x) at a point x0
     * using the symmetric difference quotient.
     * f'(x0) ≈ (f(x0 + h) - f(x0 - h)) / (2h)
     * @tparam T The arithmetic type of the function's input and output.
     * @tparam Func The type of the callable function (e.g., lambda, std::function).
     * @param f The function for which to approximate the derivative.
     * @param x0 The point at which to approximate the derivative.
     * @param h A small positive value representing the step size for the difference quotient.
     * A default of 1e-6 is chosen for common floating-point precision.
     * @return An optional containing the approximated derivative, or std::nullopt if 'h' is
     * not positive or if function evaluations at (x0 + h) or (x0 - h) result in
     * NaN or Inf.
     */
    export template <Arithmetic T, typename Func>
        std::optional<T> derivative(Func f, T x0, T h = static_cast<T>(1e-6)) {
        if (h <= 0) {
            // Step size 'h' must be positive.
            return std::nullopt;
        }

        T f_plus = f(x0 + h);
        T f_minus = f(x0 - h);

        // Check for NaN or Inf results from function evaluations
        if (std::isnan(f_plus) || std::isinf(f_plus) ||
            std::isnan(f_minus) || std::isinf(f_minus)) {
            return std::nullopt; // Derivative cannot be numerically approximated from non-finite values.
        }

        return (f_plus - f_minus) / (2 * h);
    }

    /**
     * @brief Calculates the y-value on the tangent line to function f at x0, for a given x.
     * The slope of the tangent line is numerically approximated using the derivative function.
     * Tangent line equation: y - y0 = m * (x - x0) => y = y0 + m * (x - x0)
     * @tparam T The arithmetic type.
     * @tparam Func The type of the callable function (e.g., lambda, std::function).
     * @param f The original function.
     * @param x0 The x-coordinate of the point of tangency.
     * @param x The x-coordinate for which to find the corresponding y-value on the tangent line.
     * @param h_derivative The step size 'h' to use for approximating the derivative at x0.
     * Defaults to 1e-6.
     * @return An optional containing the y-value on the tangent line, or std::nullopt if
     * the derivative at x0 cannot be approximated or if f(x0) is not a finite value.
     */
    export template <Arithmetic T, typename Func>
        std::optional<T> tangent_line(Func f, T x0, T x, T h_derivative = static_cast<T>(1e-6)) {
        // Numerically approximate the slope of the tangent line at x0
        std::optional<T> slope_opt = derivative(f, x0, h_derivative);

        // If derivative could not be approximated, we cannot find the tangent line.
        if (!slope_opt.has_value()) {
            return std::nullopt;
        }

        T slope = slope_opt.value();

        // Evaluate the function at x0 to find the y-coordinate of the tangency point (y0)
        T y0 = f(x0);

        // Check if f(x0) itself is a finite number
        if (std::isnan(y0) || std::isinf(y0)) {
            return std::nullopt; // Cannot form tangent line if point of tangency is non-finite.
        }

        // Calculate the y-value on the tangent line
        return y0 + slope * (x - x0);
    }

} // namespace calculus