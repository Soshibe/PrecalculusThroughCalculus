/*
    7. Systems of Equations and Inequalities
        Solving linear systems algebraically and graphically

        Systems of nonlinear equations

        Systems in three variables

        Matrix operations and inverses

        Determinants and Cramer's Rule
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>

import systems_equations;

int main() {
    using namespace systems;
    using T = double;
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "=== Solving Linear Systems ===\n";

    // Solve 2x2 system:
    // 3x + 4y = 25
    // 5x + 2y = 7
    Mat2<T> A2{ {{3, 4}, {5, 2}} };
    Vec2<T> b2{ 25, 7 };
    auto sol2 = solve_2x2(A2, b2);
    if (sol2) {
        std::cout << "2x2 System Solution: (" << sol2->x << ", " << sol2->y << ")\n";
    }
    else {
        std::cout << "2x2 System has no unique solution.\n";
    }

    // Solve 3x3 system:
    // x + 2y + 3z = 14
    // 4x + 5y + 6z = 32
    // 7x + 8y + 9z = 50 (note: determinant will be zero here)
    Mat3<T> A3{ {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}} };
    Vec3<T> b3{ 14, 32, 50 };
    auto sol3 = solve_3x3(A3, b3);
    if (sol3) {
        std::cout << "3x3 System Solution: (" << sol3->x << ", " << sol3->y << ", " << sol3->z << ")\n";
    }
    else {
        std::cout << "3x3 System has no unique solution.\n";
    }

    std::cout << "\n=== Matrix Inverses and Determinants ===\n";

    auto inv2 = A2.inverse();
    std::cout << "Determinant 2x2: " << A2.det() << "\n";
    if (inv2) {
        std::cout << "Inverse 2x2:\n";
        for (int i = 0; i < 2; ++i) {
            std::cout << inv2->m[i][0] << " " << inv2->m[i][1] << "\n";
        }
    }
    else {
        std::cout << "2x2 Matrix is singular.\n";
    }

    std::cout << "Determinant 3x3: " << A3.det() << "\n";
    auto inv3 = A3.inverse();
    if (inv3) {
        std::cout << "Inverse 3x3:\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << inv3->m[i][0] << " " << inv3->m[i][1] << " " << inv3->m[i][2] << "\n";
        }
    }
    else {
        std::cout << "3x3 Matrix is singular.\n";
    }

    std::cout << "\n=== Cramer's Rule ===\n";

    auto cramer2 = cramer_2x2(A2, b2);
    if (cramer2) {
        std::cout << "Cramer 2x2: (" << cramer2->x << ", " << cramer2->y << ")\n";
    }
    else {
        std::cout << "Cramer 2x2 has no unique solution.\n";
    }

    auto cramer3 = cramer_3x3(A3, b3);
    if (cramer3) {
        std::cout << "Cramer 3x3: (" << cramer3->x << ", " << cramer3->y << ", " << cramer3->z << ")\n";
    }
    else {
        std::cout << "Cramer 3x3 has no unique solution.\n";
    }

    std::cout << "\n=== Systems of Nonlinear Equations ===\n";

    // Example nonlinear system:
    // f1(x,y) = x^2 + y^2 - 25 = 0 (circle radius 5)
    // f2(x,y) = y - x = 0 (line y = x)

    auto f1 = [](T x, T y) { return x * x + y * y - 25; };
    auto f2 = [](T x, T y) { return y - x; };

    // Search grid [-6, 6] in both x and y, 50 steps each
    auto nonlinear_roots = solve_nonlinear_2var<T>(f1, f2, -6.0, 6.0, -6.0, 6.0);

    if (!nonlinear_roots.empty()) {
        for (size_t i = 0; i < nonlinear_roots.size(); ++i) {
            std::cout << "Root " << (i + 1) << ": ("
                << nonlinear_roots[i].x << ", "
                << nonlinear_roots[i].y << ")\n";
        }
    }
    else {
        std::cout << "No roots found in the specified grid.\n";
    }
    char c;
    std::cout << "Press any key to exit...\n";
    std::cin.get(c); // Wait for user input before exiting
    return 0;
}
