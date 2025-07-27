import calculus_function_mapper;

#include <iostream>
#include <string>
#include <map>
#include <functional>

std::string to_string(x_ContinuityType type) {
    switch (type) {
    case x_ContinuityType::Continuous: return "Continuous";
    case x_ContinuityType::PointDiscontinuity: return "Point Discontinuity";
    case x_ContinuityType::JumpDiscontinuity: return "Jump Discontinuity";
    case x_ContinuityType::InfiniteDiscontinuity: return "Infinite Discontinuity";
    case x_ContinuityType::RemovableDiscontinuity: return "Removable Discontinuity";
    default: return "Unknown";
    }
}

int main() {
    constexpr double xmin = -2.0;
    constexpr double xmax = 2.0;
    constexpr double ymin = -10.0;
    constexpr double ymax = 10.0;

    const int res_level = 10;

    // f2: Step function — jump discontinuity at x = 0
    auto f2 = [](double x) { return x < 0 ? -1.0 : 1.0; };

    // f3: 1 / x — infinite discontinuity at x = 0
    auto f3 = [](double x) { return 1.0 / x; };

    // f4: x^2 — continuous everywhere
    auto f4 = [](double x) { return x * x; };

    // f5: x^3 — continuous, steep slope near origin
    auto f5 = [](double x) { return x * x * x; };
    
    // f6: sin(1/x), oscillatory at x = 0
    auto f6 = [](double x) { return x == 0.0 ? 0.0 : std::sin(1.0 / x); };

    // f7: x * sin(1/x), continuous at x = 0
    auto f7 = [](double x) { return x == 0.0 ? 0.0 : x * std::sin(1.0 / x); };

    // f8: sin(1/x) / x, infinite at x = 0
    auto f8 = [](double x) { return x == 0.0 ? 0.0 : std::sin(1.0 / x) / x; };


    std::map<std::string, std::function<double(double)>> functions = {
        { "f2 (step function)", f2 },
        { "f3 (1/x)", f3 },
        { "f4 (x^2)", f4 },
        { "f5 (x^3)", f5 },
        { "f6 (sin(1/x))", f6 },
        { "f7 (x * sin(1/x))", f7 },
		{ "f8 (sin(1/x) / x)", f8 }
    };

    for (const auto& [name, func] : functions) {
        std::cout << "Analyzing: " << name << "\n";

        FunctionMapper<decltype(func)> mapper(func, xmin, xmax, ymin, ymax, res_level, false);
        mapper.map_function();
        x_ContinuityType type = mapper.check_continuity();

        std::cout << "→ Detected Continuity: " << to_string(type) << "\n\n";
    }
    std::cin.get();
    return 0;
}
