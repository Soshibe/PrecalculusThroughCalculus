export module sequences_probability;
import <vector>;
import <cmath>;
import <optional>;
import <concepts>;

export namespace sequences {

    template <typename T>
    concept Arithmetic = std::is_arithmetic_v<T>;

    // Arithmetic Sequence: aₙ = a₁ + (n - 1)d
    export template <Arithmetic T>
        T arithmetic_nth(T a1, T d, int n) {
        return a1 + (n - 1) * d;
    }

    // Arithmetic Series: Sₙ = n/2 * (2a₁ + (n-1)d)
    export template <Arithmetic T>
        T arithmetic_sum(T a1, T d, int n) {
        // Option 1 (preferable for mixed types):
        return (static_cast<T>(n) / 2.0) * (2 * a1 + (n - 1) * d); // Explicitly use 2.0 for float/double
        // Option 2 (if T is guaranteed float/double):
        return (static_cast<T>(n) / 2) * (2 * a1 + (n - 1) * d);
        // Option 3 (if T is integer and you want exact integer arithmetic result where possible)
        // This is often not how series sums are conceptualized if intermediate floats are involved.
        // It's usually S_n = n * (a1 + an) / 2
        // S_n = (n * (2 * a1 + (n - 1) * d)) / 2; (careful with overflow for large T)
    }

    // Geometric Sequence: aₙ = a₁ * r^(n - 1)
    export template <Arithmetic T>
        T geometric_nth(T a1, T r, int n) {
        return a1 * std::pow(r, n - 1);
    }

    // Geometric Series: Sₙ = a₁(1 - rⁿ) / (1 - r) for r ≠ 1
    export template <Arithmetic T>
        std::optional<T> geometric_sum(T a1, T r, int n) {
        if (r == 1) return std::nullopt;
        return a1 * (1 - std::pow(r, n)) / (1 - r);
    }

    // Infinite Geometric Series: S = a / (1 - r), |r| < 1
    export template <Arithmetic T>
        std::optional<T> geometric_sum_infinite(T a1, T r) {
        if (std::abs(r) >= 1) return std::nullopt;
        return a1 / (1 - r);
    }

    // Sigma Notation: ∑_{k=a}^{b} f(k)
    export template <Arithmetic T, typename Func>
        T sigma(Func f, int a, int b) {
        T sum = 0;
        for (int k = a; k <= b; ++k)
            sum += f(k);
        return sum;
    }

    // Factorial
    export template <Arithmetic T>
        std::optional<T> factorial(int n) {
        if (n < 0) {
            // Factorial is not defined for negative numbers.
            return std::nullopt;
        }
        if (n == 0) {
            return static_cast<T>(1);
        }

        T result = 1;

        // Special handling for integral types to check for overflow
        if constexpr (std::is_integral_v<T>) {
            for (int i = 2; i <= n; ++i) {
                // Check for potential overflow BEFORE multiplication
                // Max value of T / i < current result means next result * i will overflow
                if (std::numeric_limits<T>::max() / i < result) {
                    return std::nullopt; // Overflow detected
                }
                result *= i;
            }
        }
        else { // For floating-point types, direct multiplication
            for (int i = 2; i <= n; ++i) {
                result *= i;
            }
        }
        return result;
    }

    // Binomial Coefficient: C(n, k) = n! / (k!(n-k)!)
    export template <Arithmetic T>
        T binomial_coefficient(int n, int k) {
        // C(n, k) = C(n, n-k), so pick smaller k
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n / 2) k = n - k; // Use C(n, k) = C(n, n-k) to reduce iterations

        T res = 1;
        for (int i = 1; i <= k; ++i) {
            res *= (n - i + 1) / i;
        }
        return res;
    }

    // Binomial Theorem Expansion: (a + b)^n
    export template <Arithmetic T>
        std::vector<T> binomial_expansion(T a, T b, int n) {
        std::vector<T> terms(n + 1);
        for (int k = 0; k <= n; ++k) {
            terms[k] = binomial_coefficient<T>(n, k) * std::pow(a, n - k) * std::pow(b, k);
        }
        return terms;
    }

    export template <Arithmetic T>
        std::optional<T> permutation(int n, int k) { // Changed return type to optional
        if (k < 0 || k > n) {
            if constexpr (std::is_integral_v<T>) {
                return static_cast<T>(0); // P(n,k) = 0 for k < 0 or k > n for integral types
            }
            else {
                return static_cast<T>(0.0); // P(n,k) = 0 for k < 0 or k > n for floating point types
            }
        }
        if (k == 0) return static_cast<T>(1); // P(n,0) = 1

        T result = 1;
        for (int i = 0; i < k; ++i) {
            // This is n * (n-1) * ... * (n-k+1)
            // Check for overflow before multiplication if T is integral
            if constexpr (std::is_integral_v<T>) {
                if (std::numeric_limits<T>::max() / (n - i) < result) {
                    return std::nullopt; // Overflow
                }
            }
            result *= (n - i);
        }
        return result;
    }

    // Basic Counting Principle: combinations nCk = n! / (k!(n - k)!)
    template <typename T>
    T combination(int n, int r) {
        if (r < 0 || r > n) return T(0);
        auto num = factorial<T>(n);
        auto den1 = factorial<T>(r);
        auto den2 = factorial<T>(n - r);
        if (!num || !den1 || !den2) return T(0);
        return *num / (*den1 * *den2);
    }

    // Mathematical Induction Simulator
    export template <Arithmetic T, typename Predicate>
        bool test_induction(Predicate P, int base, int limit) {
        for (int n = base; n <= limit; ++n)
            if (!P(n)) return false;
        return true;
    }

}
