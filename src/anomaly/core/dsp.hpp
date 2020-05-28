#pragma once

#include <eigen3/Eigen/Dense>
#include "utils.hpp"

namespace anomaly::core::dsp {

namespace detail {

// TODO(sw) replace by Eigen 3.4 slices
struct Interval {
    long m_left;
    long m_size;

    Interval(long left, long size)
        : m_left {left}
        , m_size{size} { }

    template <class T>
    Interval(const T& vector)
        : m_left{0}
        , m_size{vector.size()} { }

    Interval operator+(long x) const {
        return {m_left + x, m_size};
    }

    Interval operator-(long x) const {
        return {m_left - x, m_size};
    }
};

Interval intersect(const Interval& a, const Interval& b) {
    return {std::max(a.m_left, b.m_left), std::min(a.m_left + a.m_size, b.m_left + b.m_size)};
}

template <typename T>
auto segment(T& input, Interval interval) {
    return input.segment(interval.m_left, interval.m_size);
}

template <typename Derived>
inline auto gaussianWeights(const Eigen::VectorBlock<Derived>& values, double middleValue, double delta) {
    return (values.array() - middleValue).unaryExpr(&utils::square) / (2.0 * utils::square(delta));
}
}

// TODO(sw) factor out windowing pattern
// TODO(sw) h should be in seconds, not in array-indices
template <typename DerivedPositions, typename DerivedValues>
auto bilateralFiltering(const Eigen::MatrixBase<DerivedPositions>& positions,
                        const Eigen::MatrixBase<DerivedValues>&    input,
                        double                                     delta_d,
                        double                                     delta_i,
                        Eigen::VectorXd::Index                     h) {
    using namespace detail;
    using Index = typename Eigen::VectorXd::Index;
    Eigen::VectorXd weights{2 * h + 1};
    Eigen::VectorXd output{input.size()};
    const Interval  window{-h, 2 * h + 1}; // center at 0

    for (Index center = 0; center < input.size(); ++center) {
        const auto intersection  = intersect(window, Interval(input) - center); // the previous center is still at 0
        const auto positions_seg = segment(positions, intersection + center);
        const auto input_seg     = segment(input, intersection + center);
        auto       weights_seg   = segment(weights, {0, intersection.m_size});
        weights_seg = gaussianWeights(positions_seg, positions[center], delta_d) * gaussianWeights(input_seg, input[center], delta_i);
        const double sum_of_weights = weights.sum();
        if (abs(sum_of_weights) < 0.001) { // TODO(sw) arbitrary threshold
            // no smoothing because of numeric instability
            output[center] = input[center];
            continue;
        }
        weights_seg /= sum_of_weights;
        output[center] = weights_seg.dot(input_seg);
    }
    return output;
}

} // end namespace anomaly::core::dsp
