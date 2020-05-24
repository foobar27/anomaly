#pragma once

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

} // end namespace anomaly::core::utils
