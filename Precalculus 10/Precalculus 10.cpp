/*
    10. Conic Sections
        Parabolas

        Ellipses

        Hyperbolas

        Rotated conics

        Identifying conic sections from general equations
*/
#include <iostream>;
#include <iomanip>;
import conic_sections;

int main() {
    using namespace conics;
    std::cout << std::fixed << std::setprecision(6);

    // --- Parabolas ---
    // y = x^2 - 4x + 3
    double a = 1.0, b = -4.0, c = 3.0;
    auto vertex = vertex_from_standard(a, b, c);
    std::cout << "=== Parabolas ===\n";
    std::cout << "Vertex of y = x^2 - 4x + 3: (" << vertex.first << ", " << vertex.second << ")\n";

    auto standard_coeffs = parabola_from_vertex(vertex.first, vertex.second, a);
    std::cout << "Reconstructed standard form: y = "
        << std::get<0>(standard_coeffs) << "x^2 + "
        << std::get<1>(standard_coeffs) << "x + "
        << std::get<2>(standard_coeffs) << "\n";

    auto [focus_x, focus_y, directrix_y] = parabola_focus_directrix(a, b, c);
    std::cout << "Focus: (" << focus_x << ", " << focus_y << "), Directrix: y = " << directrix_y << "\n\n";

    // --- Ellipses ---
    Ellipse<double> ellipse{ 0.0, 0.0, 5.0, 3.0, true };
    std::cout << "=== Ellipses ===\n";
    std::cout << "Ellipse centered at (" << ellipse.h << ", " << ellipse.k << "), a = "
        << ellipse.a << ", b = " << ellipse.b
        << ", major axis along " << (ellipse.horizontal_major_axis ? "x-axis" : "y-axis") << "\n";
    auto e_ellipse = eccentricity_ellipse(ellipse);
    if (e_ellipse) {
        std::cout << "Eccentricity = " << *e_ellipse << "\n\n";
    }
    else {
        std::cout << "Eccentricity calculation invalid.\n\n";
    }

    // --- Hyperbolas ---
    Hyperbola<double> hyperbola{ 0.0, 0.0, 4.0, 3.0, true };
    std::cout << "=== Hyperbolas ===\n";
    std::cout << "Hyperbola centered at (" << hyperbola.h << ", " << hyperbola.k << "), a = "
        << hyperbola.a << ", b = " << hyperbola.b
        << ", transverse axis along " << (hyperbola.horizontal_transverse_axis ? "x-axis" : "y-axis") << "\n";
    auto e_hyperbola = eccentricity_hyperbola(hyperbola);
    if (e_hyperbola) {
        std::cout << "Eccentricity = " << *e_hyperbola << "\n\n";
    }
    else {
        std::cout << "Eccentricity calculation invalid.\n\n";
    }

    // --- Rotated Conics ---
    double A = 3.0, B = 4.0, C = 2.0;
    auto angle_rad = rotation_angle(A, B, C);
    std::cout << "=== Rotated Conics ===\n";
    std::cout << "Rotation angle for conic with A=" << A << ", B=" << B << ", C=" << C << ": "
        << angle_rad * 180.0 / PI << " degrees\n\n";

    // --- Identify Conic Section ---
    GeneralConic<double> g1{ 1.0, 0.0, 1.0, 0.0, 0.0, -4.0 }; // Circle: x^2 + y^2 = 4
    GeneralConic<double> g2{ 1.0, 0.0, -1.0, 0.0, 0.0, -1.0 }; // Hyperbola: x^2 - y^2 = 1
    GeneralConic<double> g3{ 0.0, 0.0, 1.0, 0.0, -2.0, -3.0 }; // Parabola: y^2 - 2y - 3 = 0 (rewritten form)
    GeneralConic<double> g4{ 4.0, 0.0, 9.0, 0.0, 0.0, -36.0 }; // Ellipse: 4x^2 + 9y^2 = 36

    auto classify_print = [](const GeneralConic<double>& g, std::string name) {
        ConicType t = identify_conic(g);
        std::string s;
        switch (t) {
        case ConicType::Circle: s = "Circle"; break;
        case ConicType::Ellipse: s = "Ellipse"; break;
        case ConicType::Parabola: s = "Parabola"; break;
        case ConicType::Hyperbola: s = "Hyperbola"; break;
        case ConicType::Degenerate: s = "Degenerate"; break;
        default: s = "Unknown"; break;
        }
        std::cout << name << ": " << s << "\n";
        };

    std::cout << "=== Identify Conic Section ===\n";
    classify_print(g1, "g1");
    classify_print(g2, "g2");
    classify_print(g3, "g3");
    classify_print(g4, "g4");

    std::cout << "\nPress any key to exit...\n";
    std::cin.get();
    return 0;
}
