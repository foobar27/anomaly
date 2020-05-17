#include <anomaly/core/lowess.hpp>

#include "catch.hpp"
#include "matchers.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

auto smoothen(const Eigen::VectorXd& positions, const Eigen::VectorXd& input, double f, int nSteps, double delta) {
    using namespace anomaly::core::lowess;
    auto            numberOfPoints = input.size();
    LowessAlgorithm algo(numberOfPoints);
    algo.lowess(positions, input, f, nSteps, delta);
    return algo.output();
}

TEST_CASE("Loess - Sample Data") {
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

TEST_CASE("Loess - Random Data") {
    using namespace std;

    int    points_per_period = 2000;
    int    periods           = 10;
    double y_scale           = 10.0;
    int    n                 = periods * points_per_period;

    random_device r;
    mt19937       gen(r());

    normal_distribution<float> d(0.0, 0.1 * y_scale);

    auto            positions = Eigen::VectorXd::LinSpaced(n, 1000.0, 2000.0);
    Eigen::VectorXd input(n);
    for (int i = 0; i < n; ++i) {
        input[i] = y_scale * sin(2.0 * M_PI * i / points_per_period) + d(gen);
    }
    auto output = smoothen(positions, input, 1.0 / periods / 4, 0, 10.0);

    ofstream f("/tmp/sinus.tsv");
    f << "x\tinput\toutput" << endl;
    for (int i = 0; i < n; ++i) {
        f << positions[i] << "\t" << input[i] << "\t" << output[i] << endl;
    }
    // TODO(sw) verify close to sinus curve
}
