/*
    5. Analytic Trigonometry
        Proving trig identities

        Sum and difference formulas

        Double-angle and half-angle identities

        Solving trigonometric equations

*/
import trig_equations; // Import your module

#include <iostream>   // For std::cout, std::cin
#include <cmath>      // For std::sqrt, std::pow, std::sin, std::cos
#include <optional>   // For std::optional
#include <type_traits>// For std::is_arithmetic_v (used in static_assert)
#include <vector>     // For std::vector
#include <iomanip>    // For std::setprecision

int main() {
    using T = double;
    // Compile-time assertion to ensure T is an arithmetic type, as required by the module's concepts.
    static_assert(std::is_arithmetic_v<T>, "T must be an arithmetic type for trigonometric operations.");

    // Construct a domain for identity verification.
    // It's designed to avoid direct hit on common asymptotes (e.g., PI/2 for tan)
    // by using small, non-integer-multiple-of-PI/2 steps.
    std::vector<T> domain;
    // Iterate from -PI to PI, with steps of PI/24 (~7.5 degrees)
    for (T x = -PI; x <= PI; x += PI / 24.0) {
        domain.push_back(x);
    }

    std::cout << "* PROVING TRIG IDENTITIES *\n";
    std::cout << std::fixed << std::setprecision(6); // Set precision for all output

    // Define a vector of trigonometric identities to verify.
    // Each identity includes its name, a lambda for its Left-Hand Side (LHS),
    // and a lambda for its Right-Hand Side (RHS).
    std::vector<TrigIdentity<T>> identities = {
        {
            // Pythagorean Identity
            "sin^2x + cos^2x = 1",
            [](T x) -> std::optional<T> {
                return std::pow(std::sin(x), 2) + std::pow(std::cos(x), 2);
            },
            [](T x) -> std::optional<T> {
                return 1.0;
            }
        },
        {
            // Power-Reduction Identity for Sine
            "sin^2(x) = (1 - cos(2x))/2",
            [](T x) -> std::optional<T> {
                return std::pow(std::sin(x), 2);
            },
            [](T x) -> std::optional<T> {
                return (1.0 - std::cos(2 * x)) / 2.0;
            }
        },
        {
            // Sine Sum Formula (with y = PI/4)
            "sin(x + y) = sinx*cosy + cosx*siny [y = pi/4]",
            [](T x) -> std::optional<T> {
                return std::sin(x + PI / 4.0);
            },
            [](T x) -> std::optional<T> {
                return std::sin(x) * std::cos(PI / 4.0) + std::cos(x) * std::sin(PI / 4.0);
            }
        },
        {
            // Cosine Double-Angle Identity
            "cos(2x) = cos^2x - sin^2x",
            [](T x) -> std::optional<T> {
                return std::cos(2 * x);
            },
            [](T x) -> std::optional<T> {
                return std::pow(std::cos(x), 2) - std::pow(std::sin(x), 2);
            }
        },
        {
            // Tangent Double-Angle Identity
            // Using tan_robust for handling asymptotes in both LHS and RHS
            "tan(2x) = 2tanx / (1 - tan^2x)",
            [](T x) -> std::optional<T> { // LHS
                return tan_robust<T>(2 * x);
            },
            [](T x) -> std::optional<T> { // RHS
                auto t = tan_robust<T>(x); // Calculate tan(x) robustly
                if (!t.has_value()) {
                    return std::nullopt; // tan(x) is undefined, so RHS is undefined
                }
                // Check for division by zero for the RHS denominator (1 - tan^2(x))
                if (std::abs(1 - t.value() * t.value()) < EPS) {
                    return std::nullopt; // Denominator is close to zero, RHS is undefined
                }
                return (2 * t.value()) / (1 - t.value() * t.value());
            }
        }
    };

    // Verify all defined identities over the specified domain.
    verify_all_identities<T>(identities, domain);

    std::cout << "\n* SOLVING TRIGONOMETRIC EQUATIONS *\n";

    // Demonstrate symbolic solutions for basic trigonometric equations.
    // The module provides the general form of the solution.
    solve_trig_equation_symbolic<T>("sin", 0.5);          // sin(x) = 1/2
    solve_trig_equation_symbolic<T>("tan", std::sqrt(3.0)); // tan(x) = sqrt(3)
    solve_trig_equation_symbolic<T>("cos", -0.707);       // cos(x) = -0.707 (approx -sqrt(2)/2)
    solve_trig_equation_symbolic<T>("sin", 2.0);          // sin(x) = 2 (no real solution)

    // Demonstrate numeric solutions for a trigonometric equation within a specific range.
    // It finds approximate roots by stepping through the domain.
    std::cout << "\n--- Numeric solutions for cos(x) = -0.707 in [-PI, PI] ---\n";
    auto numeric_roots_cos = solve_trig_equation_numeric<T>(
        [](T x) -> std::optional<T> {
            return std::cos(x);
        },
        -0.707, -PI, PI // Target value, min x, max x
    );
    print_numeric_solutions(numeric_roots_cos); // Print the found roots

    std::cout << "\n--- Numeric solutions for tan(x) = 1.0 in [-PI, PI] ---\n";
    auto numeric_roots_tan = solve_trig_equation_numeric<T>(
        [](T x) -> std::optional<T> {
            return tan_robust<T>(x); // Use robust tan for numerical solving
        },
        1.0, -PI, PI
    );
    print_numeric_solutions(numeric_roots_tan);

    // Keep console window open until user input
    std::cout << "\nPress any key to exit...";
    char c;
    std::cin >> c;

    return 0;
}