/*
8. Sequences, Series, and Probability
    Arithmetic sequences and series
    Geometric sequences and series
    Sigma notation
    Mathematical induction
    Binomial theorem
    Basics of probability and counting
*/

import sequences_probability;
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <optional>

using std::cout;
using std::endl;
using std::optional;
using std::vector;

int main() {
    cout << std::fixed << std::setprecision(6);
    using T = double;

    // === Arithmetic Sequences and Series ===
    T a1 = 3.0, d = 2.5;
    int n = 10;

    cout << "=== Arithmetic Sequences and Series ===\n";
    cout << "10th term: " << sequences::arithmetic_nth<T>(a1, d, n) << endl;
    cout << "Sum of first 10 terms: " << sequences::arithmetic_sum<T>(a1, d, n) << "\n\n";

    // === Geometric Sequences and Series ===
    T r = 2.0;

    cout << "=== Geometric Sequences and Series ===\n";
    cout << "10th term: " << sequences::geometric_nth<T>(a1, r, n) << endl;
    auto geoSum = sequences::geometric_sum<T>(a1, r, n);
    if (geoSum) cout << "Sum of first 10 terms: " << *geoSum << endl;
    auto geoInf = sequences::geometric_sum_infinite<T>(a1, 0.5);
    if (geoInf) cout << "Infinite geometric sum (r = 0.5): " << *geoInf << "\n\n";

    // === Sigma Notation ===
    cout << "=== Sigma Notation ===\n";
    auto sumSquares = sequences::sigma<T>([](int k) { return k * k; }, 1, 5);
    cout << "Sum of k^2 from 1 to 5: " << sumSquares << "\n\n";

    // === Mathematical Induction Simulator ===
    cout << "=== Mathematical Induction ===\n";
    auto P = [](int n) -> bool {
        T lhs = sequences::sigma<T>([](int k) { return k; }, 1, n);  // sum of first n natural numbers
        T rhs = n * (n + 1) / 2.0;
        return std::abs(lhs - rhs) < 1e-6;
        };
    cout << "P(n): 1 + 2 + ... + n = n(n+1)/2 holds for n = 1 to 20? "
        << (sequences::test_induction<T>(P, 1, 20) ? "True" : "False") << "\n\n";

    // === Binomial Theorem ===
    cout << "=== Binomial Theorem Expansion ===\n";
    T a = 2.0, b = 3.0;
    int exp = 5;
    vector<T> expansion = sequences::binomial_expansion<T>(a, b, exp);
    cout << "(2 + 3)^5 expansion terms:\n";
    for (int i = 0; i <= exp; ++i)
        cout << "Term " << i << ": " << expansion[i] << endl;
    cout << "\n";

    // === Basics of Probability and Counting ===
    cout << "=== Probability and Counting ===\n";

    auto fac = sequences::factorial<T>(5);
    if (fac) cout << "5! = " << *fac << endl;

    auto perm = sequences::permutation<T>(10, 3);
    if (perm) cout << "P(10, 3) = " << *perm << endl;

    cout << "C(10, 3) = " << sequences::combination<T>(10, 3) << "\n";
    char c;
    std::cout << "Press any key to exit...\n";
    std::cin.get(c); // Wait for user input before exiting
    return 0;
}
