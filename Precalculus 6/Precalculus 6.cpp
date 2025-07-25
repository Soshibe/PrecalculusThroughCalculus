/*

6. Applications of Trigonometry
Law of Sines and Cosines

Area of triangles

Solving triangles

Trigonometric form of complex numbers

Polar coordinates and graphs
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <optional>
#include <string>
#include <utility>   // for std::pair
#include <concepts>

import trig_applications;
import trig_equations; 

int main() {
    using T = double;
    std::cout << std::fixed << std::setprecision(6);

    // --- LAW OF SINES AND COSINES ---
    std::cout << "=== Law of Sines and Cosines ===\n";

    // Law of Sines: Given side a=7, angle A=40°, find side c opposite angle C=70°
    T a = 7;
    T A = deg_to_rad(40);
    T C = deg_to_rad(70);
    auto c_opt = law_of_sines_side(a, A, C);
    if (c_opt) std::cout << "Law of Sines side c: " << c_opt.value() << "\n";

    // Law of Cosines: Given sides a=8, b=6 and angle C=60°, find side c
    a = 8; T b = 6;
    T C_cos = deg_to_rad(60);
    auto c_cos_opt = law_of_cosines_side(a, b, C_cos);
    if (c_cos_opt) std::cout << "Law of Cosines side c: " << c_cos_opt.value() << "\n";

    // --- AREA OF TRIANGLES ---
    std::cout << "\n=== Area of Triangles ===\n";

    // Area using SAS formula: sides 8 and 6 with included angle 60°
    auto area_sas = triangle_area_sas(a, b, C_cos);
    if (area_sas) std::cout << "Area (SAS): " << area_sas.value() << "\n";

    // Area using Heron's formula: sides 7, 8, 9
    T c_heron = 9;
    auto area_heron = triangle_area_heron(7, 8, 9);
    if (area_heron) std::cout << "Area (Heron's): " << area_heron.value() << "\n";

    // --- SOLVING TRIANGLES ---
    std::cout << "\n=== Solving Triangles ===\n";

    // Solve SSS triangle with sides 7, 8, 9
    auto tri_sss = solve_triangle_sss(7, 8, 9);
    if (tri_sss) {
        const auto& t = tri_sss.value();
        std::cout << "SSS Triangle angles (deg): A=" << rad_to_deg(t.A)
            << ", B=" << rad_to_deg(t.B)
            << ", C=" << rad_to_deg(t.C) << "\n";
    }

    // Solve SAS triangle: sides a=7, b=8 and included angle C=50°
    auto tri_sas = solve_triangle_sas(7, deg_to_rad(50), 8);
    if (tri_sas) {
        const auto& t = tri_sas.value();
        std::cout << "SAS Triangle sides: a=" << t.a << ", b=" << t.b << ", c=" << t.c << "\n";
        std::cout << "SAS Triangle angles (deg): A=" << rad_to_deg(t.A)
            << ", B=" << rad_to_deg(t.B)
            << ", C=" << rad_to_deg(t.C) << "\n";
    }

    // Solve ASA triangle: angles 40°, 60° and side a=7
    auto tri_asa = solve_triangle_asa_aas(deg_to_rad(40), deg_to_rad(60), 7);
    if (tri_asa) {
        const auto& t = tri_asa.value();
        std::cout << "ASA Triangle sides: a=" << t.a << ", b=" << t.b << ", c=" << t.c << "\n";
        std::cout << "ASA Triangle angles (deg): A=" << rad_to_deg(t.A)
            << ", B=" << rad_to_deg(t.B)
            << ", C=" << rad_to_deg(t.C) << "\n";
    }

    // Solve SSA ambiguous case: a=7, b=10, A=30°
    auto ssa_solutions = solve_triangle_ssa(7, 10, deg_to_rad(30));
    std::cout << "SSA ambiguous case solutions found: " << ssa_solutions.size() << "\n";
    for (size_t i = 0; i < ssa_solutions.size(); ++i) {
        const auto& t = ssa_solutions[i];
        std::cout << "Solution " << (i + 1) << ": sides a=" << t.a << ", b=" << t.b << ", c=" << t.c << "\n";
        std::cout << "Angles (deg): A=" << rad_to_deg(t.A)
            << ", B=" << rad_to_deg(t.B)
            << ", C=" << rad_to_deg(t.C) << "\n";
    }

    // --- TRIGONOMETRIC FORM OF COMPLEX NUMBERS ---
    std::cout << "\n=== Trigonometric Form of Complex Numbers ===\n";
    T x = 3, y = 4;
    auto polar = to_polar(x, y);
    std::cout << "Cartesian (" << x << ", " << y << ") -> Polar (r=" << polar.r << ", theta="
        << rad_to_deg(polar.theta) << " degrees)\n";

    auto cart = to_cartesian(polar);
    std::cout << "Polar -> Cartesian (" << cart.first << ", " << cart.second << ")\n";

    // Multiplying complex numbers in polar form (3+4i) * (1+i)
    auto p1 = polar;
    auto p2 = to_polar(1, 1);
    auto p_mul = multiply_polar(p1, p2);
    auto cart_mul = to_cartesian(p_mul);
    std::cout << "Multiplication result (polar): r=" << p_mul.r << ", theta=" << rad_to_deg(p_mul.theta) << " degrees\n";
    std::cout << "Multiplication result (cartesian): (" << cart_mul.first << ", " << cart_mul.second << ")\n";

    // --- POLAR COORDINATES AND GRAPHS ---
    std::cout << "\n=== Polar Coordinates ===\n";
    T r = 5;
    T theta = deg_to_rad(45);
    auto cart_polar = polar_to_cartesian(r, theta);
    std::cout << "Polar (r=5, theta=45°) -> Cartesian (" << cart_polar.first << ", " << cart_polar.second << ")\n";

    T x_p = 3, y_p = 3;
    auto polar_cart = cartesian_to_polar(x_p, y_p);
    std::cout << "Cartesian (3,3) -> Polar (r=" << polar_cart.first << ", theta=" << rad_to_deg(polar_cart.second) << " degrees)\n";

    std::cout << "\nPress any key to exit...";
    char c;
    std::cin >> c;

    return 0;
}
