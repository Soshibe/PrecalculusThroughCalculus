
#include <functional>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <utility>
#include <optional>
#include <set>
#include <memory>
#include <limits>
#include <concepts>
#include <type_traits>
#include <iomanip>
export module calculus_function_mapper;

export enum class x_ContinuityType {
    Continuous,
    PointDiscontinuity,
    JumpDiscontinuity,
    InfiniteDiscontinuity,
    RemovableDiscontinuity,
    Unknown
};

export template <typename F>
class FunctionPredictor {
public:
    FunctionPredictor(F func, double x_min, double x_max, int resolution)
        : f(func), x_min(x_min), x_max(x_max), resolution(resolution) {

        if (resolution <= 0) resolution = 1;
        dx = (x_max - x_min) / resolution;

        if (std::abs(dx) < std::numeric_limits<double>::epsilon())
            dx = 1.0;

        x_last = x_max;
        f_last = f(x_last);

        double x_prev = x_last - dx;
        slope_last = std::abs(dx) < std::numeric_limits<double>::epsilon() ?
            0.0 : (f_last - f(x_prev)) / dx;
    }

    double predictNext() const {
        return f_last + slope_last * dx;
    }

    static double get_local_slope(F func, double x1, double x2) {
        if (std::abs(x2 - x1) < std::numeric_limits<double>::epsilon())
            return std::numeric_limits<double>::infinity();
        return (func(x2) - func(x1)) / (x2 - x1);
    }

private:
    F f;
    double x_min, x_max, dx;
    int resolution;
    double x_last, f_last, slope_last;
};

export template<typename Func, typename T_Domain = double, typename T_Output = double>
class FunctionMapper {
public:
    FunctionMapper(Func func,
        T_Domain x_min, T_Domain x_max,
        T_Domain y_min, T_Domain y_max,
        int resolution_level = 9,
        bool quiet_mode = false);

    void map_function();
    x_ContinuityType check_continuity() const;
    const std::map<int, std::optional<int>>& get_mapped_grid() const { return mapped_grid_; }

private:
    Func func_;
    T_Domain x_min_, x_max_, y_min_, y_max_;
    int resolution_level_;
    bool quiet_mode_;

    int grid_width_;
    int grid_height_;

    std::map<int, std::optional<int>> mapped_grid_;
    std::set<int> infinite_discontinuity_x_coords_;

    int scale_y_to_grid(T_Output y_val) const;

    void print_message(const std::string& msg) const {
        if (!quiet_mode_) std::cout << msg << '\n';
    }
};

template<typename Func, typename T_Domain, typename T_Output>
FunctionMapper<Func, T_Domain, T_Output>::FunctionMapper(Func func,
    T_Domain x_min, T_Domain x_max,
    T_Domain y_min, T_Domain y_max,
    int resolution_level,
    bool quiet_mode)
    : func_(std::move(func)),
    x_min_(x_min), x_max_(x_max),
    y_min_(y_min), y_max_(y_max),
    resolution_level_(resolution_level),
    quiet_mode_(quiet_mode) {

    if (x_min_ >= x_max_ || y_min_ >= y_max_) {
        print_message("Error: Invalid domain.");
        grid_width_ = grid_height_ = 0;
        return;
    }

    const T_Domain dx = x_max_ - x_min_;
    const T_Domain dy = y_max_ - y_min_;
    long long multiplier = 0;

    if (resolution_level_ >= -1 && resolution_level_ <= 61) {
        multiplier = 1LL << (resolution_level_ + 1);
    }
    else {
        print_message("Warning: Invalid resolution level. Using fallback level 9.");
        resolution_level_ = 9;
        multiplier = 1LL << (resolution_level_ + 1);
    }

    grid_width_ = static_cast<int>(dx * multiplier) + 1;
    grid_height_ = static_cast<int>(dy * multiplier) + 1;
    if (grid_width_ <= 0) grid_width_ = 1;
    if (grid_height_ <= 0) grid_height_ = 1;

    print_message("Initialized mapper grid: " + std::to_string(grid_width_) + "x" + std::to_string(grid_height_));
}

template<typename Func, typename T_Domain, typename T_Output>
int FunctionMapper<Func, T_Domain, T_Output>::scale_y_to_grid(T_Output y_val) const {
    if (y_max_ == y_min_) return 0;
    T_Output norm = (y_val - y_min_) / (y_max_ - y_min_);
    return static_cast<int>(std::round(norm * (grid_height_ - 1)));
}

template<typename Func, typename T_Domain, typename T_Output>
void FunctionMapper<Func, T_Domain, T_Output>::map_function() {
    mapped_grid_.clear();
    infinite_discontinuity_x_coords_.clear();

    if (grid_width_ == 0 || grid_height_ == 0) {
        print_message("Cannot map function: Grid dimensions zero.");
        return;
    }

    T_Domain x_step = (x_max_ - x_min_) / static_cast<T_Domain>(grid_width_ - 1);

    for (int gx = 0; gx < grid_width_; ++gx) {
        T_Domain x = x_min_ + gx * x_step;
        T_Output y = func_(x);

        if (std::isinf(y) || std::isnan(y)) {
            mapped_grid_[gx] = std::nullopt;
            if (std::isinf(y)) infinite_discontinuity_x_coords_.insert(gx);
        }
        else {
            int gy = scale_y_to_grid(y);
            if (gy >= 0 && gy < grid_height_) mapped_grid_[gx] = gy;
            else mapped_grid_[gx] = std::nullopt;
        }
    }

    print_message("Mapped " + std::to_string(mapped_grid_.size()) + " columns.");
}

template<typename Func, typename T_Domain, typename T_Output>
x_ContinuityType FunctionMapper<Func, T_Domain, T_Output>::check_continuity() const {
    if (mapped_grid_.empty()) {
        print_message("Continuity check failed: grid not mapped.");
        return x_ContinuityType::Unknown;
    }

    x_ContinuityType result = x_ContinuityType::Continuous;
    bool jump = false, inf = false, removable = false;

    for (int gx = 0; gx < grid_width_; ++gx) {
        auto cur_it = mapped_grid_.find(gx);
        auto cur = (cur_it != mapped_grid_.end()) ? cur_it->second : std::nullopt;

        auto prev_it = (gx > 0) ? mapped_grid_.find(gx - 1) : mapped_grid_.end();
        auto prev = (prev_it != mapped_grid_.end()) ? prev_it->second : std::nullopt;

        auto next_it = (gx < grid_width_ - 1) ? mapped_grid_.find(gx + 1) : mapped_grid_.end();
        auto next = (next_it != mapped_grid_.end()) ? next_it->second : std::nullopt;

        if (infinite_discontinuity_x_coords_.contains(gx)) {
            inf = true;
            result = x_ContinuityType::InfiniteDiscontinuity;
            continue;
        }

        if (!cur.has_value() && prev.has_value() && next.has_value()
            && std::abs(prev.value() - next.value()) <= 1) {
            removable = true;
            if (result == x_ContinuityType::Continuous)
                result = x_ContinuityType::RemovableDiscontinuity;
            continue;
        }

        bool jump_candidate = false;
        if (!cur.has_value() && (prev.has_value() || next.has_value())) {
            jump_candidate = true;
        }
        else if (cur.has_value() && prev.has_value() &&
            std::abs(cur.value() - prev.value()) > 1) {
            jump_candidate = true;
        }

        if (jump_candidate && gx >= 2) {
            T_Domain x_step = (x_max_ - x_min_) / static_cast<T_Domain>(grid_width_ - 1);
            T_Domain x0 = x_min_ + (gx - 2) * x_step;
            T_Domain x1 = x_min_ + (gx - 1) * x_step;
            T_Domain x2 = x_min_ + gx * x_step;

            FunctionPredictor<Func> pred(func_, x0, x1, 1);
            T_Output predicted_y = pred.predictNext();
            T_Output actual_y = func_(x2);

            T_Output tol = std::max((y_max_ - y_min_) * 0.01, std::numeric_limits<T_Output>::epsilon() * 100);
            if (std::abs(actual_y - predicted_y) > tol) {
                jump = true;
                if (result == x_ContinuityType::Continuous || result == x_ContinuityType::RemovableDiscontinuity)
                    result = x_ContinuityType::JumpDiscontinuity;
            }
        }
    }

    // Multi-resolution validation
    if (jump || inf) {
        print_message("--- Discontinuity detected, verifying with lower resolutions ---");
        bool confirmed = true;

        for (int r = resolution_level_ - 1; r >= -1; --r) {
            FunctionMapper<Func, T_Domain, T_Output> sub(func_, x_min_, x_max_, y_min_, y_max_, r, true);
            sub.map_function();
            x_ContinuityType res = sub.check_continuity();

            if (res != x_ContinuityType::JumpDiscontinuity &&
                res != x_ContinuityType::InfiniteDiscontinuity) {
                confirmed = false;
                break;
            }
        }

        if (confirmed) {
            return inf ? x_ContinuityType::InfiniteDiscontinuity : x_ContinuityType::JumpDiscontinuity;
        }
        else {
            print_message("--- Not confirmed across all resolutions. Reclassifying...");
            return removable ? x_ContinuityType::RemovableDiscontinuity : x_ContinuityType::Continuous;
        }
    }

    return result;
}
