#include <anomaly/core/stl.hpp>

#include "catch.hpp"
#include "matchers.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace anomaly::core::stl;

static auto decompose(StlConfiguration config, const Eigen::VectorXd& input) {
    StlAlgorithm algorithm(config, input.size());
    algorithm.stl(input);
    return algorithm;
}

TEST_CASE("STL - Sinus with noise") {
    using namespace std;

    int    points_per_period = 2000;
    int    periods           = 10;
    double y_scale           = 10.0;
    int    n                 = periods * points_per_period;

    random_device r;
    mt19937       gen(r());

    normal_distribution<double> d(0.0, 0.1 * y_scale);

    Eigen::VectorXd input(n);
    for (int i = 0; i < n; ++i) {
        input[i] = y_scale * sin(2.0 * M_PI * i / points_per_period) + d(gen);
        // TODO(sw) add artificial trend, e.g. parabola
    }
    StlConfiguration config{points_per_period,
                            // seasonal smoother
                            {7, 1, 0.0, 1, 1},
                            // trend smoother
                            {points_per_period * 2, 1, 0.0, 1, 1},
                            // lowpass smoother
                            {2001, 1, 0.0, 1, 1},
                            // #iterations inner loop
                            1,
                            // #iterations outer loop
                            1};

    auto decomposition = decompose(config, input);

    ofstream f("/tmp/sinus-decomposed.tsv");
    f << "x\tinput\tseason\ttrend" << endl;
    for (int i = 0; i < n; ++i) {
        f << i << "\t" << input[i] << "\t" << decomposition.season()[i] << "\t" << decomposition.trend()[i] << endl;
    }
    // TODO(sw) verify close to sinus curve
}
