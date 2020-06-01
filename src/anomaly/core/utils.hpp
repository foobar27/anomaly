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
concept StaticOrDynamicSize = requires {
    Value >= -1;
};

template <Eigen::Index Value>
concept StaticSize = requires {
    Value >= 0;
};

template <Eigen::Index Value>
concept DynamicSize = requires {
    Value == -1;
};

template <Eigen::Index Value> requires StaticOrDynamicSize<Value>
struct StaticOrDynamicSizeContainer {
    static constexpr Eigen::Index staticValue = Value;
    static constexpr bool isStatic = true;
    static_assert(Value >= 0);

    constexpr StaticOrDynamicSizeContainer() { }

    constexpr auto operator()() const {
        return Value;
    }

//    template<Eigen::Index OtherValue>
//    constexpr StaticOrDynamicSizeContainer<Value + OtherValue> operator+(constexpr OtherValue) const requires StaticSize<OtherValue>{
//        return {};
//    }
};

template <>
struct StaticOrDynamicSizeContainer<Eigen::Dynamic> {
    static constexpr Eigen::Index staticValue = Eigen::Dynamic;
    static constexpr bool isStatic = false;

    constexpr StaticOrDynamicSizeContainer(Eigen::Index value)
        : m_value(value) { }

    auto operator()() const {
        return m_value;
    }

private:
    Eigen::Index m_value;
};

} // end namespace anomaly::core::utils
