#pragma once

#include <eigen3/Eigen/Dense>

namespace anomaly::core::utils {

constexpr inline double square(double x) {
    return x * x;
}

constexpr inline double biSquare(double x) {
    return square(1.0 - square(x));
}

constexpr inline double cube(double x) {
    return x * x * x;
}

constexpr inline double triCube(double x) {
    return cube(1.0 - cube(x));
}

template <Eigen::Index Value>
struct StaticOrDynamicSize {
    static constexpr bool isStatic    = true;
    static constexpr Eigen::Index    staticValue = Value;
    static_assert(Value >= 0);

    StaticOrDynamicSize() { }

    constexpr Eigen::Index operator()() const {
        return Value;
    }

};

template <>
struct StaticOrDynamicSize<Eigen::Dynamic> {
    static constexpr bool isStatic    = false;
    static constexpr Eigen::Index    staticValue = Eigen::Dynamic;

    StaticOrDynamicSize(Eigen::Index value)
        : m_value(value) { }

    Eigen::Index operator()() const {
        return m_value;
    }

private:
    Eigen::Index m_value;
};

} // end namespace anomaly::core::utils
