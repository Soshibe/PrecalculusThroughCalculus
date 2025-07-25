export module systems_equations;

import <vector>;
import <optional>;
import <cmath>;
import <functional>; // Required for std::function
import <concepts>;   // Required for concepts

export namespace systems {

    // ADDED 'inline' KEYWORD HERE
    export inline constexpr double EPS_LINEAR = 1e-9;
    export inline constexpr double EPS_NONLINEAR_ROOT = 1e-3;

    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    template <Arithmetic T>
    struct Vec2 {
        T x, y;
    };

    template <Arithmetic T>
    struct Vec3 {
        T x, y, z;
    };

    template <Arithmetic T>
    struct Mat2 {
        T m[2][2];

        T det() const {
            return m[0][0] * m[1][1] - m[0][1] * m[1][0];
        }

        std::optional<Mat2<T>> inverse() const {
            T d = det();
            if (std::abs(d) < EPS_LINEAR) return std::nullopt;
            Mat2<T> inv;
            inv.m[0][0] = m[1][1] / d;
            inv.m[0][1] = -m[0][1] / d;
            inv.m[1][0] = -m[1][0] / d;
            inv.m[1][1] = m[0][0] / d;
            return inv;
        }
    };

    template <Arithmetic T>
    struct Mat3 {
        T m[3][3];

        T det() const {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }

        std::optional<Mat3<T>> inverse() const {
            T d = det();
            if (std::abs(d) < EPS_LINEAR) return std::nullopt;
            Mat3<T> inv;

            inv.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
            inv.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d;
            inv.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;

            inv.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d;
            inv.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
            inv.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d;

            inv.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
            inv.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d;
            inv.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;

            return inv;
        }
    };

    export template <Arithmetic T>
        std::optional<Vec2<T>> solve_2x2(const Mat2<T>& A, const Vec2<T>& b) {
        auto A_inv = A.inverse();
        if (!A_inv) return std::nullopt;
        const auto& m = A_inv.value().m;
        return Vec2<T>{
            m[0][0] * b.x + m[0][1] * b.y,
                m[1][0] * b.x + m[1][1] * b.y
        };
    }

    export template <Arithmetic T>
        std::optional<Vec3<T>> solve_3x3(const Mat3<T>& A, const Vec3<T>& b) {
        auto A_inv = A.inverse();
        if (!A_inv) return std::nullopt;
        const auto& m = A_inv.value().m;
        return Vec3<T>{
            m[0][0] * b.x + m[0][1] * b.y + m[0][2] * b.z,
                m[1][0] * b.x + m[1][1] * b.y + m[1][2] * b.z,
                m[2][0] * b.x + m[2][1] * b.y + m[2][2] * b.z
        };
    }

    export template <Arithmetic T>
        std::optional<Vec2<T>> cramer_2x2(const Mat2<T>& A, const Vec2<T>& b) {
        T det_A = A.det();
        if (std::abs(det_A) < EPS_LINEAR) return std::nullopt;

        Mat2<T> A1 = A;
        A1.m[0][0] = b.x; A1.m[1][0] = b.y;
        T det_A1 = A1.det();

        Mat2<T> A2 = A;
        A2.m[0][1] = b.x; A2.m[1][1] = b.y;
        T det_A2 = A2.det();

        return Vec2<T>{ det_A1 / det_A, det_A2 / det_A };
    }

    export template <Arithmetic T>
        std::optional<Vec3<T>> cramer_3x3(const Mat3<T>& A, const Vec3<T>& b) {
        T det_A = A.det();
        if (std::abs(det_A) < EPS_LINEAR) return std::nullopt;

        Mat3<T> A_x = A;
        A_x.m[0][0] = b.x; A_x.m[1][0] = b.y; A_x.m[2][0] = b.z;
        T det_Ax = A_x.det();

        Mat3<T> A_y = A;
        A_y.m[0][1] = b.x; A_y.m[1][1] = b.y; A_y.m[2][1] = b.z;
        T det_Ay = A_y.det();

        Mat3<T> A_z = A;
        A_z.m[0][2] = b.x; A_z.m[1][2] = b.y; A_z.m[2][2] = b.z;
        T det_Az = A_z.det();

        return Vec3<T>{ det_Ax / det_A, det_Ay / det_A, det_Az / det_A };
    }

    export template <Arithmetic T>
        std::optional<Vec2<T>> newton_refine_2var(
            std::function<T(T, T)> f1,
            std::function<T(T, T)> f2,
            Vec2<T> start,
            int max_iters = 20,
            T tol = 1e-9,
            T h = 1e-6)
    {
        T x = start.x;
        T y = start.y;

        for (int iter = 0; iter < max_iters; ++iter) {
            T F1 = f1(x, y);
            T F2 = f2(x, y);

            if (std::abs(F1) < tol && std::abs(F2) < tol)
                return Vec2<T>{ x, y };

            // Central difference for partial derivatives
            T dF1dx = (f1(x + h, y) - f1(x - h, y)) / (2 * h);
            T dF1dy = (f1(x, y + h) - f1(x, y - h)) / (2 * h);
            T dF2dx = (f2(x + h, y) - f2(x - h, y)) / (2 * h);
            T dF2dy = (f2(x, y + h) - f2(x, y - h)) / (2 * h);

            T det_jacobian = dF1dx * dF2dy - dF1dy * dF2dx;
            if (std::abs(det_jacobian) < 1e-12) return std::nullopt; // Jacobian is singular or near-singular

            // Solving for delta_x and delta_y using Cramer's rule on the linear system J * delta_x = -F
            // dx_step = (dF2dy * F1 - dF1dy * F2) / det_jacobian;  // This was effectively -dx
            // dy_step = (dF1dx * F2 - dF2dx * F1) / det_jacobian;  // This was effectively -dy
            // Correct update for Newton's method: x_new = x_old - delta_x
            // delta_x = J_inv * F
            // The expressions you had were for -delta_x, so it correctly subtracts directly.
            T dx_step = (dF2dy * F1 - dF1dy * F2) / det_jacobian;
            T dy_step = (dF1dx * F2 - dF2dx * F1) / det_jacobian;

            x -= dx_step;
            y -= dy_step;

            if (std::sqrt(dx_step * dx_step + dy_step * dy_step) < tol)
                return Vec2<T>{ x, y };
        }
        return std::nullopt; // No convergence within max_iters
    }

    export template <Arithmetic T>
        std::vector<Vec2<T>> solve_nonlinear_2var(
            std::function<T(T, T)> f1,
            std::function<T(T, T)> f2,
            T x_start, T x_end,
            T y_start, T y_end)
    {
        const T dx_span = std::abs(x_end - x_start);
        const T dy_span = std::abs(y_end - y_start);

        // 🔁 Dynamic resolution based on range size
        constexpr T resolution = 0.01; // ~100 points per unit
        const int x_steps = std::max(10, static_cast<int>(dx_span / resolution));
        const int y_steps = std::max(10, static_cast<int>(dy_span / resolution));

        T dx = dx_span / x_steps;
        T dy = dy_span / y_steps;

        std::vector<Vec2<T>> candidates;

        // 🧲 Loose detection pass
        constexpr T detection_eps = 0.2;

        for (int i = 0; i < x_steps; ++i) {
            for (int j = 0; j < y_steps; ++j) {
                T x0 = x_start + i * dx;
                T y0 = y_start + j * dy;
                T x1 = x0 + dx;
                T y1 = y0 + dy;

                // Sample corners and center
                T f1_a = f1(x0, y0), f1_b = f1(x1, y0), f1_c = f1(x0, y1), f1_d = f1(x1, y1);
                T f2_a = f2(x0, y0), f2_b = f2(x1, y0), f2_c = f2(x0, y1), f2_d = f2(x1, y1);
                T f1_m = f1((x0 + x1) / 2, (y0 + y1) / 2);
                T f2_m = f2((x0 + x1) / 2, (y0 + y1) / 2);

                bool f1_loose = std::min({ std::abs(f1_a), std::abs(f1_b), std::abs(f1_c), std::abs(f1_d), std::abs(f1_m) }) < detection_eps;
                bool f2_loose = std::min({ std::abs(f2_a), std::abs(f2_b), std::abs(f2_c), std::abs(f2_d), std::abs(f2_m) }) < detection_eps;

                if (f1_loose && f2_loose) {
                    Vec2<T> guess = { (x0 + x1) / 2, (y0 + y1) / 2 };

                    bool duplicate = false;
                    for (const auto& pt : candidates) {
                        if (std::hypot(pt.x - guess.x, pt.y - guess.y) < dx * 0.5) {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                        candidates.push_back(guess);
                }
            }
        }

        // ✨ Newton refinement for final polish
        constexpr T tol = 1e-7;
        std::vector<Vec2<T>> refined;

        for (const auto& guess : candidates) {
            auto root = newton_refine_2var(f1, f2, guess);
            if (root) {
                bool is_duplicate = false;
                for (const auto& r : refined) {
                    if (std::hypot(r.x - root->x, r.y - root->y) < tol) {
                        is_duplicate = true;
                        break;
                    }
                }
                if (!is_duplicate)
                    refined.push_back(*root);
            }
        }

        return refined;
    }

} // namespace systems