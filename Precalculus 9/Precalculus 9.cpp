/*
    9. Vectors and Parametric Equations
        Vector operations

        Dot product

        Applications of vectors

        Parametric equations and their graphs

        Eliminating the parameter
*/
import vectors_parametric;
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace vectors;

int main() {
    std::cout << std::fixed << std::setprecision(6);

    Vec2<double> a{ 3.0, 4.0 };
    Vec2<double> b{ 1.0, 2.0 };
    double scalar = 2.5;

    // --- Vector Operations ---
    std::cout << "=== Vector Operations ===\n";
    auto a_plus_b = add(a, b);
    auto a_scaled = scale(a, scalar);
    std::cout << "a + b = (" << a_plus_b.first << ", " << a_plus_b.second << ")\n";
    std::cout << "a * " << scalar << " = (" << a_scaled.first << ", " << a_scaled.second << ")\n";
    std::cout << "|a| = " << magnitude(a) << "\n\n";

    // --- Dot Product ---
    std::cout << "=== Dot Product and Angle ===\n";
    auto dot_ab = dot(a, b);
    std::cout << "a · b = " << dot_ab << "\n";
    auto angle = angle_between(a, b);
    if (angle)
        std::cout << "Angle between a and b = " << (*angle * 180.0 / PI) << " degrees\n\n";
    else
        std::cout << "Angle undefined (zero vector)\n\n";

    // --- Vector Projection ---
    std::cout << "=== Projection of a onto b ===\n";
    auto proj = projection(a, b);
    if (proj)
        std::cout << "proj_b(a) = (" << proj->first << ", " << proj->second << ")\n\n";
    else
        std::cout << "Cannot project onto zero vector\n\n";

    // --- Parametric Equations ---
    std::cout << "=== Parametric Equations ===\n";
    Vec2<double> pos0{ 1.0, 2.0 };
    Vec2<double> velocity{ 3.0, -1.0 };
    double t = 4.0;
    auto pos_t = parametric_position(pos0, velocity, t);
    std::cout << "Position at t = " << t << ": (" << pos_t.first << ", " << pos_t.second << ")\n\n";

    // --- Eliminate Parameter ---
    std::cout << "=== Eliminate Parameter ===\n";
    auto result = eliminate_parameter(pos0, velocity);
    if (result) {
        auto [m, b] = *result;
        std::cout << "y = " << m << "x + " << b << "\n";
    }
    else {
        std::cout << "Cannot eliminate parameter (velocity.x = 0)\n";
    }

    std::cout << "\nPress any key to exit...\n";
    std::cin.get();
}
