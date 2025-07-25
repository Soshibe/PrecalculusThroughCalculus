import limits_intro;

#include <iostream>
#include <iomanip>

using namespace calc::limits;

int main() {
    using x_Real = double;

    auto x_f1 = [](x_Real x) { return (x * x - 1) / (x - 1); };
    x_Real x_a1 = 1.0;

    x_LimitApproximator<x_Real, decltype(x_f1)> x_approx1;

    std::cout << "=== Limit Approximation ===\n";
    auto x_L = x_approx1.x_left_limit(x_f1, x_a1);
    auto x_R = x_approx1.x_right_limit(x_f1, x_a1);
    auto x_T = x_approx1.x_two_sided_limit(x_f1, x_a1);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Left limit at x = " << x_a1 << ": " << (x_L ? *x_L : NAN) << '\n';
    std::cout << "Right limit at x = " << x_a1 << ": " << (x_R ? *x_R : NAN) << '\n';
    std::cout << "Two-sided limit at x = " << x_a1 << ": " << (x_T ? *x_T : NAN) << "\n\n";

    auto x_f2 = [](x_Real x) { return (x == 2.0) ? 5.0 : x * x; };
    x_LimitApproximator<x_Real, decltype(x_f2)> x_approx2;
    std::cout << "=== Discontinuity Classification ===\n";
    std::cout << "f(x) = x^2, but f(2) = 5\n";
    std::cout << "Discontinuity type at x = 2: "
        << x_approx2.x_discontinuity_type(x_f2, 2.0) << "\n\n";

    auto x_f3 = [](x_Real x) { return (x < 0) ? -1.0 : 1.0; };
    x_LimitApproximator<x_Real, decltype(x_f3)> x_approx3;
    std::cout << "f(x) = -1 for x < 0, 1 for x >= 0\n";
    std::cout << "Discontinuity type at x = 0: "
        << x_approx3.x_discontinuity_type(x_f3, 0.0) << "\n\n";

    auto x_f4 = [](x_Real x) { return 1.0 / x; };
    x_LimitApproximator<x_Real, decltype(x_f4)> x_approx4;
    std::cout << "f(x) = 1/x\n";
    std::cout << "Discontinuity type at x = 0: "
        << x_approx4.x_discontinuity_type(x_f4, 0.0) << "\n\n";

    std::cout << "=== Sample Table Near x = 1 for f(x) = (x^2 - 1)/(x - 1) ===\n";
    auto x_table = x_approx1.x_table_near(x_f1, x_a1);
    for (const auto& x_s : x_table) {
        std::cout << "x = " << x_s.x_x << ", f(x) = " << (x_s.x_fx ? *x_s.x_fx : NAN) << '\n';
    }
    std::cin.get();
    return 0;
}
