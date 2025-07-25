/*
    11. Limits and an Introduction to Calculus (optional in some texts)
        Understanding limits intuitively

        Limits graphically and numerically

        Introduction to the concept of a derivative

        Tangent lines
*/

#include <iostream>
#include <iomanip>
#include <cmath>
import limits_calculus;

int main() {
    using namespace calculus;
    std::cout << std::fixed << std::setprecision(6);

    // Function for limit and derivative tests: f(x) = (x^2 - 1) / (x - 1)
    auto f = [](double x) -> double {
        if (x == 1.0) return 2.0; // Define f(1) by simplifying (x^2 - 1)/(x-1) = x+1 except at x=1
        return (x * x - 1) / (x - 1);
        };

    double point = 1.0;

    // --- Limit Approximation ---
    auto limit_res = approximate_limit(f, point);
    std::cout << "=== Limit Approximation ===\n";
    if (limit_res.left_limit)
        std::cout << "Left limit at x = " << point << ": " << *limit_res.left_limit << '\n';
    else
        std::cout << "Left limit at x = " << point << ": undefined\n";

    if (limit_res.right_limit)
        std::cout << "Right limit at x = " << point << ": " << *limit_res.right_limit << '\n';
    else
        std::cout << "Right limit at x = " << point << ": undefined\n";

    if (limit_res.has_two_sided_limit())
        std::cout << "Two-sided limit at x = " << point << ": " << *limit_res.get_two_sided_limit() << "\n\n";
    else
        std::cout << "Two-sided limit at x = " << point << " does not exist or does not agree.\n\n";

    // --- Derivative Approximation ---
    auto deriv = derivative(f, point);
    std::cout << "=== Derivative Approximation ===\n";
    if (deriv)
        std::cout << "Approximate derivative of f at x = " << point << ": " << *deriv << "\n\n";
    else
        std::cout << "Derivative at x = " << point << " could not be approximated.\n\n";

    // --- Tangent Line ---
    std::cout << "=== Tangent Line ===\n";
    double x_val = 1.5;
    auto tangent_y = tangent_line(f, point, x_val);
    if (tangent_y) {
        std::cout << "Value of tangent line to f at x = " << point << " evaluated at x = " << x_val << ": " << *tangent_y << "\n";
    }
    else {
        std::cout << "Could not compute tangent line at x = " << point << "\n";
    }
    char c;
	std::cin >> c; // Wait for user input before exiting
    return 0;
}
