#include <vector>
#include <concepts>
#include <limits>
#include <string>
#include <optional>
#include <ranges>
#include <algorithm>
export module limits_intro;
export namespace calc::limits {

    namespace x_detail {
        template <typename T>
        constexpr bool x_is_nan(T x_val) {
            return x_val != x_val;
        }

        template <typename T>
        constexpr bool x_is_inf(T x_val) {
            return x_val == std::numeric_limits<T>::infinity() ||
                x_val == -std::numeric_limits<T>::infinity();
        }

        template <typename T>
        constexpr T x_abs(T x_val) {
            return (x_val < static_cast<T>(0)) ? -x_val : x_val;
        }

        template <typename T>
        constexpr T x_clamp(T x_val, T x_min_val, T x_max_val) {
            return (x_val < x_min_val) ? x_min_val : (x_val > x_max_val) ? x_max_val : x_val;
        }
    }

    export template <typename T>
        concept x_FloatingPoint = std::floating_point<T>;

    export template <x_FloatingPoint T>
        struct x_Sample {
        T x_x;
        std::optional<T> x_fx;
    };

    export template <x_FloatingPoint T, typename x_Func>
        requires std::invocable<x_Func, T>&&
    std::convertible_to<std::invoke_result_t<x_Func, T>, T>
        struct x_LimitApproximator {
        static constexpr T x_delta_scale = static_cast<T>(1e-6);
        static constexpr T x_delta_min = static_cast<T>(1e-12);
        static constexpr T x_delta_max = static_cast<T>(1e-2);

        static constexpr T x_epsilon_scale = static_cast<T>(1e-5);
        static constexpr T x_epsilon_min = static_cast<T>(1e-12);
        static constexpr T x_epsilon_max = static_cast<T>(1e-2);

        T x_adaptive_delta(T x_a) const {
            T x_base_delta = x_detail::x_abs(x_a) * x_delta_scale;
            return x_detail::x_clamp(x_base_delta, x_delta_min, x_delta_max);
        }

        T x_adaptive_epsilon(x_Func const& x_f, T x_a) const {
            auto x_val_opt = x_safe_eval(x_f, x_a);
            T x_val = x_val_opt ? x_detail::x_abs(*x_val_opt) : static_cast<T>(1);
            T x_base_epsilon = x_val * x_epsilon_scale;
            return x_detail::x_clamp(x_base_epsilon, x_epsilon_min, x_epsilon_max);
        }

        std::optional<T> x_safe_eval(x_Func const& x_f, T x_val) const {
            T x_result = x_f(x_val);
            if (x_detail::x_is_nan(x_result) || x_detail::x_is_inf(x_result)) return std::nullopt;
            return x_result;
        }

        std::optional<T> x_left_limit(x_Func const& x_f, T x_a) const {
            T x_delta = x_adaptive_delta(x_a);
            return x_safe_eval(x_f, x_a - x_delta);
        }

        std::optional<T> x_right_limit(x_Func const& x_f, T x_a) const {
            T x_delta = x_adaptive_delta(x_a);
            return x_safe_eval(x_f, x_a + x_delta);
        }

        std::optional<T> x_two_sided_limit(x_Func const& x_f, T x_a) const {
            T x_delta = x_adaptive_delta(x_a);
            T x_epsilon = x_adaptive_epsilon(x_f, x_a);
            auto x_L = x_safe_eval(x_f, x_a - x_delta);
            auto x_R = x_safe_eval(x_f, x_a + x_delta);
            if (x_L && x_R && x_detail::x_abs(*x_L - *x_R) <= x_epsilon) {
                return (*x_L + *x_R) / static_cast<T>(2);
            }
            return std::nullopt;
        }

        std::vector<x_Sample<T>> x_table_near(x_Func const& x_f, T x_a,
            T x_range = static_cast<T>(0.01),
            int x_steps = 5) const {
            T x_step_size = x_range / x_steps;
            auto x_indices = std::views::iota(1, x_steps + 1);
            std::vector<x_Sample<T>> x_samples;

            std::ranges::for_each(x_indices | std::views::reverse, [&](int x_i) {
                T x_val = x_a - x_i * x_step_size;
                x_samples.emplace_back(x_Sample<T>{x_val, x_safe_eval(x_f, x_val)});
                });

            std::ranges::for_each(x_indices, [&](int x_i) {
                T x_val = x_a + x_i * x_step_size;
                x_samples.emplace_back(x_Sample<T>{x_val, x_safe_eval(x_f, x_val)});
                });

            return x_samples;
        }

        std::string x_discontinuity_type(x_Func const& x_f, T x_a) const {
            T x_delta = x_adaptive_delta(x_a);
            T x_epsilon = x_adaptive_epsilon(x_f, x_a);
            auto x_L = x_safe_eval(x_f, x_a - x_delta);
            auto x_R = x_safe_eval(x_f, x_a + x_delta);
            auto x_fx = x_safe_eval(x_f, x_a);

            if (!x_L || !x_R || !x_fx) return "Undefined/Oscillating discontinuity";
            if (x_detail::x_is_inf(*x_L) || x_detail::x_is_inf(*x_R) || x_detail::x_is_inf(*x_fx)) {
                return "Infinite discontinuity";
            }
            if (x_detail::x_abs(*x_L - *x_R) > x_epsilon) {
                return "Jump discontinuity";
            }
            if (x_detail::x_abs(*x_fx - *x_L) > x_epsilon) {
                return "Removable discontinuity";
            }
            return "Continuous";
        }
    };

} // namespace calc::limits