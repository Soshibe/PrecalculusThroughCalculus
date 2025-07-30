import <concepts>;
import <vector>;
import <functional>;
import <cmath>;
import <optional>;
import <limits>;
import <sstream>;
import <iomanip>;
import <algorithm>;

export module limits_formal_full_stateful_domain;

export template <typename Real>
concept Arithmetic = std::is_arithmetic_v<Real>;

export enum class LimitCheckResult {
    HoldsForAllSampledEpsilons,
    FailsForSomeEpsilons,
    LimitDoesNotExist,
    Undefined
};

export template <Arithmetic Real>
class EpsilonDeltaLimit {
private:
    Real a;                      // Approach point (can be ±∞)
    Real L;                      // Candidate limit
    std::function<std::optional<Real>(Real)> f;

    std::vector<Real> epsilon_samples;
    std::vector<Real> delta_candidates;

    std::optional<Real> domain_min;
    std::optional<Real> domain_max;

    LimitCheckResult cached_result = LimitCheckResult::Undefined;
    std::vector<std::optional<Real>> cached_deltas; // per epsilon

    bool oscillating_behavior_cached = false;
    bool left_fail_cached = false;
    bool right_fail_cached = false;

    // --- New cached slope/intercept for infinite limits ---
    mutable std::optional<std::pair<Real, Real>> cached_linear_asymptote; // slope, intercept

    bool is_infinite_a() const {
        return std::isinf(a);
    }

    // Helper: sample function values between x_start and x_end with num_samples
    std::vector<std::optional<Real>> sample_function(Real x_start, Real x_end, size_t num_samples) const {
        std::vector<std::optional<Real>> samples;
        if (num_samples == 0) return samples;
        Real step = (x_end - x_start) / static_cast<Real>(num_samples - 1);
        for (size_t i = 0; i < num_samples; ++i) {
            Real x = x_start + i * step;
            if ((domain_min && x < *domain_min) || (domain_max && x > *domain_max)) {
                samples.push_back(std::nullopt);
            }
            else {
                samples.push_back(f(x));
            }
        }
        return samples;
    }

    // Estimate linear asymptote y = mx + b on interval [x_start, x_end]
    std::optional<std::pair<Real, Real>> estimate_linear_asymptote(Real x_start, Real x_end, size_t samples = 10) const {
        auto vals = sample_function(x_start, x_end, samples);
        // Check for missing values
        for (auto& v : vals) if (!v.has_value()) return std::nullopt;

        Real sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
        for (size_t i = 0; i < vals.size(); ++i) {
            Real x = x_start + (x_end - x_start) * static_cast<Real>(i) / static_cast<Real>(vals.size() - 1);
            Real y = vals[i].value();
            sum_x += x;
            sum_y += y;
            sum_xx += x * x;
            sum_xy += x * y;
        }
        Real n = static_cast<Real>(vals.size());
        Real denom = n * sum_xx - sum_x * sum_x;
        if (denom == 0) return std::nullopt;
        Real m = (n * sum_xy - sum_x * sum_y) / denom;
        Real b = (sum_y - m * sum_x) / n;
        return std::make_pair(m, b);
    }

    // Check epsilon-delta for finite a
    bool epsilon_delta_condition_finite(Real epsilon, Real delta) const {
        Real left_x = a - delta;
        Real right_x = a + delta;
        if ((domain_min && left_x < *domain_min) || (domain_max && left_x > *domain_max))
            return false;
        if ((domain_min && right_x < *domain_min) || (domain_max && right_x > *domain_max))
            return false;
        for (Real dx : {left_x, right_x}) {
            if (dx == a) continue;
            auto fx = f(dx);
            if (!fx.has_value() || std::abs(fx.value() - L) >= epsilon) {
                return false;
            }
        }
        return true;
    }

    // Check epsilon-delta for infinite a (limit at ±∞)
    bool epsilon_delta_condition_infinite(Real epsilon, Real M) const {
        std::vector<Real> test_points;
        if (a > 0) {
            test_points = { M, M + std::abs(M) * 0.1 };
        }
        else {
            test_points = { M, M - std::abs(M) * 0.1 };
        }
        for (Real x : test_points) {
            if ((domain_min && x < *domain_min) || (domain_max && x > *domain_max))
                return false;
            auto fx = f(x);
            if (!fx.has_value() || std::abs(fx.value() - L) >= epsilon)
                return false;
        }
        return true;
    }

    // Check linear asymptote near infinity for zero slope & intercept ~ L
    bool check_linear_asymptote_near_infinity(Real epsilon) const {
        if (!cached_linear_asymptote.has_value()) {
            Real large_M = delta_candidates.empty() ? static_cast<Real>(1e6) : delta_candidates.back();

            Real x_start = (a > 0) ? (large_M) : (-large_M * 10);
            Real x_end = (a > 0) ? (large_M * 10) : (-large_M);

            cached_linear_asymptote = estimate_linear_asymptote(x_start, x_end, 50);
            if (!cached_linear_asymptote.has_value()) return false;
        }

        Real m = cached_linear_asymptote->first;
        Real b = cached_linear_asymptote->second;

        if (std::abs(m) <= epsilon) {
            return (std::abs(b - L) <= epsilon);
        }
        return false;
    }

    // Find position where second derivative of average slope <= epsilon, toward infinity
    std::optional<Real> find_position_for_slope_decay(Real start, Real step, size_t max_samples, Real epsilon) const {
        if (max_samples < 5) return std::nullopt; // Need enough samples for second derivative

        std::vector<Real> xs(max_samples);
        for (size_t i = 0; i < max_samples; ++i) {
            xs[i] = (a > 0) ? (start + i * step) : (start - i * step);
            if ((domain_min && xs[i] < *domain_min) || (domain_max && xs[i] > *domain_max))
                return std::nullopt; // Out of domain
        }

        // Compute average slopes m_i between xs[i] and xs[i+1]
        std::vector<Real> m;
        for (size_t i = 0; i + 1 < xs.size(); ++i) {
            auto f1 = f(xs[i]);
            auto f2 = f(xs[i + 1]);
            if (!f1.has_value() || !f2.has_value())
                return std::nullopt;
            m.push_back((f2.value() - f1.value()) / (xs[i + 1] - xs[i]));
        }

        if (m.size() < 3) return std::nullopt;

        // First derivative of slope dm_i = (m[i+1] - m[i]) / (xs[i+2] - xs[i])
        std::vector<Real> dm;
        for (size_t i = 0; i + 2 < xs.size(); ++i) {
            dm.push_back((m[i + 1] - m[i]) / (xs[i + 2] - xs[i]));
        }

        if (dm.size() < 2) return std::nullopt;

        // Second derivative d2m_i = (dm[i+1] - dm[i]) / (xs[i+3] - xs[i])
        for (size_t i = 0; i + 3 < xs.size(); ++i) {
            Real d2m = (dm[i + 1] - dm[i]) / (xs[i + 3] - xs[i]);
            if (std::abs(d2m) <= epsilon) {
                return xs[i + 1]; // position where slope decay <= epsilon
            }
        }

        return std::nullopt;
    }

    void compute_verification() {
        cached_deltas.clear();
        bool all_eps_pass = true;

        for (Real epsilon : epsilon_samples) {
            bool found_delta = false;

            if (is_infinite_a()) {
                // First try linear asymptote zero slope case
                if (check_linear_asymptote_near_infinity(epsilon)) {
                    for (auto d : delta_candidates) {
                        cached_deltas.push_back(d);
                        found_delta = true;
                    }
                    break;
                }
                else {
                    // Try to find position where slope decay matches epsilon
                    std::optional<Real> slope_decay_pos = find_position_for_slope_decay(
                        delta_candidates.empty() ? static_cast<Real>(1e3) : delta_candidates.front(),
                        0.1 * (delta_candidates.empty() ? static_cast<Real>(1e3) : delta_candidates.front()),
                        20,
                        epsilon
                    );

                    Real start_pos = slope_decay_pos.value_or(
                        delta_candidates.empty() ? static_cast<Real>(1e3) : delta_candidates.front()
                    );

                    // Now test delta candidates starting from this position
                    for (Real delta : delta_candidates) {
                        if (delta < start_pos) continue;
                        if (!epsilon_delta_condition_infinite(epsilon, delta)) continue;
                        cached_deltas.push_back(delta);
                        found_delta = true;
                        break;
                    }
                }
            }
            else {
                // Finite a usual check
                for (Real delta : delta_candidates) {
                    if (delta <= 0) continue;
                    if ((domain_min && (a - delta < *domain_min)) ||
                        (domain_max && (a + delta > *domain_max))) {
                        continue;
                    }
                    if (!epsilon_delta_condition_finite(epsilon, delta))
                        continue;
                    cached_deltas.push_back(delta);
                    found_delta = true;
                    break;
                }
            }

            if (!found_delta) {
                cached_deltas.push_back(std::nullopt);
                cached_result = LimitCheckResult::FailsForSomeEpsilons;
                return;
            }
        }

        cached_result = all_eps_pass ? LimitCheckResult::HoldsForAllSampledEpsilons : LimitCheckResult::Undefined;
    }

    void compute_oscillating_behavior() {
        Real min_epsilon = *std::min_element(epsilon_samples.begin(), epsilon_samples.end());
        Real base_r = std::min(min_epsilon, static_cast<Real>(1e-6));
        if (is_infinite_a()) {
            Real test_M = delta_candidates.empty() ? static_cast<Real>(1e6) : delta_candidates.front();
            Real x1 = (a > 0) ? test_M : -test_M;
            Real x2 = (a > 0) ? (test_M + base_r) : (-test_M - base_r);
            if ((domain_min && x1 < *domain_min) || (domain_max && x1 > *domain_max)) {
                oscillating_behavior_cached = false;
                return;
            }
            if ((domain_min && x2 < *domain_min) || (domain_max && x2 > *domain_max)) {
                oscillating_behavior_cached = false;
                return;
            }
            auto f1 = f(x1);
            auto f2 = f(x2);
            if (!f1.has_value() || !f2.has_value()) {
                oscillating_behavior_cached = false;
                return;
            }
            Real delta_y = std::abs(f1.value() - f2.value());
            oscillating_behavior_cached = (delta_y > min_epsilon);
        }
        else {
            Real left_limit = domain_min.value_or(a - base_r);
            Real right_limit = domain_max.value_or(a + base_r);
            Real r = base_r;
            if (a - r < left_limit) r = a - left_limit;
            if (a + r > right_limit) r = right_limit - a;
            if (r <= 0) {
                oscillating_behavior_cached = false;
                return;
            }
            Real x1 = a - r;
            Real x2 = a + r;
            if ((domain_min && x1 < *domain_min) || (domain_max && x1 > *domain_max)) {
                oscillating_behavior_cached = false;
                return;
            }
            if ((domain_min && x2 < *domain_min) || (domain_max && x2 > *domain_max)) {
                oscillating_behavior_cached = false;
                return;
            }
            auto f1 = f(x1);
            auto f2 = f(x2);
            if (!f1.has_value() || !f2.has_value()) {
                oscillating_behavior_cached = false;
                return;
            }
            Real delta_y = std::abs(f1.value() - f2.value());
            oscillating_behavior_cached = (delta_y > min_epsilon);
        }
    }

    void compute_left_fail() {
        if (is_infinite_a()) {
            left_fail_cached = false;
            return;
        }
        for (Real epsilon : epsilon_samples) {
            bool left_fails = true;
            for (Real delta : delta_candidates) {
                Real x_left = a - delta;
                if ((domain_min && x_left < *domain_min) || (domain_max && x_left > *domain_max)) continue;
                auto fx = f(x_left);
                if (fx.has_value() && std::abs(fx.value() - L) < epsilon) {
                    left_fails = false;
                    break;
                }
            }
            if (left_fails) {
                left_fail_cached = true;
                return;
            }
        }
        left_fail_cached = false;
    }

    void compute_right_fail() {
        if (is_infinite_a()) {
            right_fail_cached = false;
            return;
        }
        for (Real epsilon : epsilon_samples) {
            bool right_fails = true;
            for (Real delta : delta_candidates) {
                Real x_right = a + delta;
                if ((domain_min && x_right < *domain_min) || (domain_max && x_right > *domain_max)) continue;
                auto fx = f(x_right);
                if (fx.has_value() && std::abs(fx.value() - L) < epsilon) {
                    right_fails = false;
                    break;
                }
            }
            if (right_fails) {
                right_fail_cached = true;
                return;
            }
        }
        right_fail_cached = false;
    }

public:
    EpsilonDeltaLimit(
        Real a_,
        Real L_,
        std::function<std::optional<Real>(Real)> f_,
        std::optional<Real> domain_min_ = std::nullopt,
        std::optional<Real> domain_max_ = std::nullopt,
        std::vector<Real> epsilon_samples_ = { 1.0, 0.1, 0.01, 0.001 },
        std::vector<Real> delta_candidates_ = { 0.1, 0.01, 0.001, 0.0001 }
    )
        : a(a_), L(L_), f(std::move(f_)),
        domain_min(domain_min_), domain_max(domain_max_),
        epsilon_samples(std::move(epsilon_samples_)),
        delta_candidates(std::move(delta_candidates_))
    {
        compute_verification();
        compute_left_fail();
        compute_right_fail();
        compute_oscillating_behavior();
    }

    LimitCheckResult result() const {
        return cached_result;
    }

    std::optional<Real> find_delta_for_epsilon(Real epsilon) const {
        for (size_t i = 0; i < epsilon_samples.size(); ++i) {
            if (epsilon_samples[i] == epsilon && i < cached_deltas.size()) {
                return cached_deltas[i];
            }
        }
        return std::nullopt;
    }

    std::string proof_report() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < epsilon_samples.size(); ++i) {
            Real epsilon = epsilon_samples[i];
            if (i < cached_deltas.size() && cached_deltas[i].has_value()) {
                if (is_infinite_a()) {
                    oss << "For epsilon = " << epsilon << ", M = " << cached_deltas[i].value()
                        << ": |f(x) - L| < epsilon whenever x ";
                    if (a > 0) oss << ">= M.\n";
                    else oss << "<= M.\n";
                }
                else {
                    oss << "For epsilon = " << epsilon << ", delta = " << cached_deltas[i].value()
                        << ": |f(x) - L| < epsilon whenever 0 < |x - a| < delta.\n";
                }
            }
            else {
                oss << "There exists epsilon = " << epsilon << " for which no delta satisfies the condition.\n";
            }
        }
        return oss.str();
    }

    bool fails_from_left_only() const {
        return left_fail_cached;
    }

    bool fails_from_right_only() const {
        return right_fail_cached;
    }

    bool has_oscillating_behavior() const {
        return oscillating_behavior_cached;
    }
};
