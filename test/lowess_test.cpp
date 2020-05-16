#include <anomaly/core/lowess.hpp>

#include "catch.hpp"
#include "matchers.hpp"

auto smoothen(const Eigen::VectorXd& positions, const Eigen::VectorXd& input, double f, int nSteps, double delta) {
    using namespace anomaly::core::lowess;
    auto            numberOfPoints = input.size();
    LowessAlgorithm algo(numberOfPoints);
    algo.lowess(positions, input, f, nSteps, delta);
    return algo.output();
}

TEST_CASE("Dummy test") {
    int numberOfPoints = 20;

    Eigen::VectorXd positions(numberOfPoints);
    positions << 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 10, 12, 14, 50;

    Eigen::VectorXd input(numberOfPoints);
    input << 18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8, 5, 19;

    SECTION("parameter set 1") {
        double          f      = 0.25;
        int             nSteps = 0;
        double          delta  = 0.0;
        Eigen::VectorXd expected(numberOfPoints);
        expected << 13.659, 11.145, 8.701, 9.722, 10.000, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300,
            13.000, 6.440, 5.596, 5.456, 18.998;

        auto actual = smoothen(positions, input, f, nSteps, delta);
        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
    }

    SECTION("parameter set 2") {
        double          f      = 0.25;
        int             nSteps = 0;
        double          delta  = 3.0;
        Eigen::VectorXd expected(numberOfPoints);
        expected << 13.659, 12.347, 11.034, 9.722, 10.511, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300,
            13.000, 6.440, 5.596, 5.456, 18.998;

        auto actual = smoothen(positions, input, f, nSteps, delta);
        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
    }

    SECTION("parameter set 3") {
        double          f      = 0.25;
        int             nSteps = 2;
        double          delta  = 0.0;
        Eigen::VectorXd expected(numberOfPoints);
        expected << 14.811, 12.115, 8.984, 9.676, 10.000, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346,
            13.000, 6.734, 5.744, 5.415, 18.998;

        auto actual = smoothen(positions, input, f, nSteps, delta);
        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
    }
}
