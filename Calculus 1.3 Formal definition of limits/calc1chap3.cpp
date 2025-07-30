// main.cpp
// Test harness for EpsilonDeltaLimit with limits at finite points and infinity

#include <iostream>
#include <optional>
#include <cmath>
#include <functional>
import limits_formal_full_stateful_domain;

using namespace std;

const char* to_string(LimitCheckResult r) {
    using R = LimitCheckResult;
    switch (r) {
    case R::HoldsForAllSampledEpsilons: return "Holds for all sampled epsilons";
    case R::FailsForSomeEpsilons: return "Fails for some epsilons";
    case R::LimitDoesNotExist: return "Limit does not exist";
    case R::Undefined: return "Undefined";
    default: return "Unknown";
    }
}

int main() {
    using Real = double;
    using EDL = EpsilonDeltaLimit<Real>;

    // Test 1: finite limit f(x) = x^2 at x=2, limit=4
    auto f1 = [](Real x) -> optional<Real> { return x * x; };
    EDL limit1(2.0, 4.0, f1);
    cout << "Test 1: f(x)=x^2 at x=2\n";
    cout << "Result: " << to_string(limit1.result()) << "\n" << limit1.proof_report() << "\n";

    // Test 2: removable discontinuity f(x)=sin(x)/x at x=0, limit=1
    auto f2 = [](Real x) -> optional<Real> {
        if (x == 0) return nullopt;
        return sin(x) / x;
        };
    EDL limit2(0.0, 1.0, f2, -1.0, 1.0);
    cout << "Test 2: f(x)=sin(x)/x at x=0\n";
    cout << "Result: " << to_string(limit2.result()) << "\n" << limit2.proof_report() << "\n";

    // Test 3: infinite discontinuity f(x)=1/x at x=0, limit does not exist
    auto f3 = [](Real x) -> optional<Real> {
        if (x == 0) return nullopt;
        return 1.0 / x;
        };
    EDL limit3(0.0, 0.0, f3, -1.0, 1.0);
    cout << "Test 3: f(x)=1/x at x=0\n";
    cout << "Result: " << to_string(limit3.result()) << "\n" << limit3.proof_report() << "\n";

    // Test 4: jump discontinuity step function at x=0
    auto f4 = [](Real x) -> optional<Real> {
        return (x < 0) ? 0.0 : 1.0;
        };
    EDL limit4(0.0, 0.0, f4, -1.0, 1.0);
    cout << "Test 4: step function at x=0\n";
    cout << "Result: " << to_string(limit4.result()) << "\n" << limit4.proof_report() << "\n";

    // Test 5: oscillating discontinuity f(x)=sin(1/x) at x=0
    auto f5 = [](Real x) -> optional<Real> {
        if (x == 0) return nullopt;
        return sin(1.0 / x);
        };
    EDL limit5(0.0, 0.0, f5, -0.1, 0.1);
    cout << "Test 5: f(x)=sin(1/x) at x=0\n";
    cout << "Result: " << to_string(limit5.result()) << "\n" << limit5.proof_report() << "\n";

    // --- Limits at infinity ---

    // Test 6: limit as x → +∞ of f(x) = 1/x, limit = 0
    auto f_inf_pos = [](Real x) -> optional<Real> {
        if (x == 0) return nullopt;
        return 1.0 / x;
        };
    EDL limit_inf_pos(
        numeric_limits<Real>::infinity(),
        0.0,
        f_inf_pos,
        1.0,           // domain_min (start testing at x=1)
        nullopt,
        { 1.0, 0.1, 0.01 },
        { 10, 100, 1000, 10000 }
    );
    cout << "Test 6: f(x)=1/x as x→+∞\n";
    cout << "Result: " << to_string(limit_inf_pos.result()) << "\n" << limit_inf_pos.proof_report() << "\n";

    // Test 7: limit as x → -∞ of f(x) = e^x, limit = 0
    auto f_inf_neg = [](Real x) -> optional<Real> {
        return exp(x);
        };
    EDL limit_inf_neg(
        -numeric_limits<Real>::infinity(),
        0.0,
        f_inf_neg,
        nullopt,
        0.0,                     // domain_max (up to 0)
        { 1.0, 0.1, 0.01 },
        { -10, -100, -1000, -10000 }
    );
    cout << "Test 7: f(x)=e^x as x->negative infinity\n";
    cout << "Result: " << to_string(limit_inf_neg.result()) << "\n" << limit_inf_neg.proof_report() << "\n";

    // --- Logarithmic function tests ---

    // Test 8: limit as x → +∞ of f(x) = ln(x), test finite L=1000 (expect fail)
    auto f_ln_inf = [](Real x) -> optional<Real> {
        if (x <= 0) return nullopt;
        return log(x);
        };
    EDL limit_ln_inf(
        numeric_limits<Real>::infinity(),
        1000.0,               // arbitrary finite limit to test failure
        f_ln_inf,
        1.0,                  // domain_min = 1 (log defined for x>0)
        nullopt,
        { 1.0, 0.1, 0.01 },
        { 10, 100, 1000, 10000 }
    );
    cout << "Test 8: f(x)=ln(x) as x→+∞ with L=1000 (expect fail)\n";
    cout << "Result: " << to_string(limit_ln_inf.result()) << "\n" << limit_ln_inf.proof_report() << "\n";

    // Test 9: limit as x → 0^+ of f(x) = ln(x), test finite L=-10 (expect fail)
    auto f_ln_zero_plus = [](Real x) -> optional<Real> {
        if (x <= 0) return nullopt;
        return log(x);
        };
    EDL limit_ln_zero_plus(
        0.0,
        -10.0,              // arbitrary finite limit to test failure
        f_ln_zero_plus,
        0.0,                // domain_min = 0, right side limit
        1.0,
        { 1.0, 0.1, 0.01 },
        { 0.1, 0.01, 0.001, 0.0001 }
    );
    cout << "Test 9: f(x)=ln(x) as x→0^+ with L=-10 (expect fail)\n";
    cout << "Result: " << to_string(limit_ln_zero_plus.result()) << "\n" << limit_ln_zero_plus.proof_report() << "\n";

    // Test 10: limit as x → 1 of f(x) = log base 2 of x, limit = 0
    auto f_log2 = [](Real x) -> optional<Real> {
        if (x <= 0) return nullopt;
        return log(x) / log(2.0);
        };
    EDL limit_log2_1(
        1.0,
        0.0,
        f_log2,
        0.1,
        2.0,
        { 1.0, 0.1, 0.01 },
        { 0.1, 0.01, 0.001 }
    );
    cout << "Test 10: f(x)=log_2(x) at x=1\n";
    cout << "Result: " << to_string(limit_log2_1.result()) << "\n" << limit_log2_1.proof_report() << "\n";
    // Test for f(x)=x*sin(x) as x->+infinity, checking for finite limit L=0 (expect fail)
    auto f_x_sinx = [](Real x) -> optional<Real> {
        return x * sin(x);
        };
    EDL limit_x_sinx_inf(
        numeric_limits<Real>::infinity(),
        0.0, // Candidate limit
        f_x_sinx,
        1.0, // domain_min (start testing at x=1)
        nullopt,
        { 1.0, 0.1, 0.01 }, // Sampled epsilons
        { 10, 100, 1000, 10000 } // Sampled M values
    );
    cout << "Test: f(x)=x*sin(x) as x->+inf with L=0 (expect fail)\n";
    cout << "Result: " << to_string(limit_x_sinx_inf.result()) << "\n" << limit_x_sinx_inf.proof_report() << "\n";
    cin.get();
    return 0;
}
