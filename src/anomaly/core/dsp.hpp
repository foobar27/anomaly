#pragma once

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <eigen3/Eigen/Dense>
#include "utils.hpp"

namespace anomaly::core::dsp {

namespace detail {

// TODO(sw) replace by Eigen 3.4 slices
struct Interval {
    using Index = Eigen::Index;
    Index m_left;
    Index m_size;

    constexpr Interval(Index left, Index size)
        : m_left{left}
        , m_size{size} { }

    template <class T>
    Interval(const T& vector)
        : m_left{0}
        , m_size{vector.size()} { }

    constexpr Interval operator+(Index x) const {
        return {m_left + x, m_size};
    }

    constexpr Interval operator-(Index x) const {
        return {m_left - x, m_size};
    }
};

constexpr Interval intersect(const Interval& a, const Interval& b) {
    auto m_left  = std::max(a.m_left, b.m_left); // inclusive
    auto m_right = std::min(a.m_left + a.m_size, b.m_left + b.m_size); // exclusive
    return {m_left, m_right - m_left};
}

template <typename T>
auto segment(T& input, Interval interval) {
    return input.segment(interval.m_left, interval.m_size);
}

}

template <typename Scalar, int WindowSize = Eigen::Dynamic>
struct MovingAverageFilter {
    using Index          = Eigen::Index;
    using WindowSizeType = utils::StaticOrDynamicSize<WindowSize>;

    template <typename = std::enable_if<WindowSizeType::isStatic>>
    MovingAverageFilter() { }

    template <typename = std::enable_if<!WindowSizeType::isStatic>>
    MovingAverageFilter(Index radius)
        : m_windowSize{radius} { }

    // If input has size n, output has size n-len+1.
    template<typename Input, typename Output>
    void operator()(const Input& input, Output& output) {
        using namespace std;
        using namespace boost::accumulators;
        accumulator_set<Scalar, stats<tag::rolling_mean>> m_acc{tag::rolling_mean::window_size = m_windowSize()};
        for_each(input.data(), input.data() + m_windowSize(), m_acc);

        // get the first average
        output[0] = rolling_mean(m_acc);

        auto new_n = input.size() - m_windowSize() + 1;
        if (new_n > 1) {
            auto k = m_windowSize() - 1;
            for (Eigen::Index j = 1; j < new_n; ++j) {
                k++;
                m_acc(input[k]);
                output[j] = rolling_mean(m_acc);
            }
        }
    }

private:
    const WindowSizeType m_windowSize;
};

template <int Radius>
struct WindowOperation {
    using Index      = Eigen::Index;
    using RadiusType = utils::StaticOrDynamicSize<Radius>;

    template <typename = std::enable_if<RadiusType::isStatic>>
    WindowOperation() { }

    template <typename = std::enable_if<!RadiusType::isStatic>>
    explicit WindowOperation(Index radius)
        : m_radius{radius} { }

    template <typename Self, typename DerivedPositions, typename DerivedValues>
    auto operator()(Self& self, const Eigen::MatrixBase<DerivedPositions>& positions, const Eigen::MatrixBase<DerivedValues>& input) const {
        using namespace detail;
        Eigen::VectorXd output{input.size()};
        const Interval  window{-m_radius(), 2 * m_radius() + 1};
        for (Index center = 0; center < input.size(); ++center) {
            const auto intersection  = intersect(window, Interval(input) - center);
            const auto positions_seg = segment(positions, intersection + center);
            const auto input_seg     = segment(input, intersection + center);
            output[center]           = self.processWindow(intersection, positions_seg, input_seg);
        }
        return output;
    }

    RadiusType m_radius;
};

template <typename Scalar, Eigen::Index Radius = Eigen::Dynamic>
struct BilateralFilter {
    using Index               = Eigen::Index;
    using WindowOperationType = WindowOperation<Radius>;
    using RadiusType          = typename WindowOperationType::RadiusType;
    using WeightsType         = Eigen::Matrix<Scalar, RadiusType::isStatic ? 2 * RadiusType::staticValue + 1 : Eigen::Dynamic, 1>;

    template <typename = std::enable_if<RadiusType::isStatic>>
    BilateralFilter() { }

    template <typename = std::enable_if<!RadiusType::isStatic>>
    explicit BilateralFilter(int radius)
        : m_windowOperation(radius)
        , m_weights(2 * radius + 1) { }

    void setDeltaD(Scalar value) {
        m_deltaD = value;
    }

    void setDeltaI(Scalar value) {
        m_deltaI = value;
    }

    template <typename DerivedPositions, typename DerivedValues>
    auto operator()(const Eigen::MatrixBase<DerivedPositions>& positions, const Eigen::MatrixBase<DerivedValues>& input) {
        return m_windowOperation(*this, positions, input);
    }

private:
    template <typename DerivedPositions, typename DerivedValues>
    auto processWindow(const detail::Interval                     window,
                       const Eigen::MatrixBase<DerivedPositions>& positions,
                       const Eigen::MatrixBase<DerivedValues>&    input) {
        using namespace detail;
        auto window_center = -window.m_left;
        auto weights       = segment(m_weights, {0, window.m_size});
        weights = gaussianWeights(positions, positions[window_center], m_deltaD) * gaussianWeights(input, input[window_center], m_deltaI);
        const auto sum_of_weights = weights.sum();
        if (abs(sum_of_weights) < 0.001) { // TODO(sw) arbitrary threshold
            // no smoothing because of numeric instability
            return input[window_center];
        }
        weights /= sum_of_weights;
        return weights.dot(input);
    }

    template <typename Derived>
    inline auto gaussianWeights(const Eigen::MatrixBase<Derived>& values, Scalar middleValue, Scalar delta) {
        return (values.array() - middleValue).unaryExpr(&utils::square) / (2.0 * utils::square(delta));
    }

    const WindowOperation<Radius> m_windowOperation;
    WeightsType                   m_weights;
    Scalar                        m_deltaD{1.0};
    Scalar                        m_deltaI{1.0};

    friend struct WindowOperation<Radius>;
};

} // end namespace anomaly::core::dsp
