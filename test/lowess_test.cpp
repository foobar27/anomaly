//#include <anomaly/core/lowess.hpp>

//#include "catch.hpp"
//#include "matchers.hpp"
//#include <cmath>
//#include <fstream>
//#include <iostream>

// static auto smoothen(const Eigen::VectorXd& positions, const Eigen::VectorXd& input, double f, int nSteps, double delta) {
//    using namespace anomaly::core::lowess;
//    auto                number_of_points = input.size();
//    LowessAlgorithm     algo(number_of_points);
//    LowessConfiguration config{0, 0, f, nSteps, delta};
//    algo.lowess(config, positions, input);
//    return algo.output();
//}

// TEST_CASE("Loess - Sample Data") {
//    int number_of_points = 20;

//    Eigen::VectorXd positions(number_of_points);
//    positions << 1, 2, 3, 4, 5, 6, 6.01, 6.02, 6.03, 6.04, 6.05, 6.06, 6.07, 6.08, 6.09, 8, 10, 12, 14, 50;

//    Eigen::VectorXd input(number_of_points);
//    input << 18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8, 5, 19;

//    SECTION("parameter set 1") {
//        double          f       = 0.25;
//        int             n_steps = 0;
//        double          delta   = 0.0;
//        Eigen::VectorXd expected(number_of_points);
//        expected << 13.6588, 11.1446, 8.70117, 9.72204, 9.99823, 9.93375, 10.1563, 11.2863, 7, 7.29466, 11.7095, 17, 16.8506, 13.9803,
//            13.1805, 12.9996, 6.37724, 5.5777, 5.5004, 18.9982;

//        auto actual = smoothen(positions, input, f, n_steps, delta);
//        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
//    }

//    SECTION("parameter set 2") {
//        double          f       = 0.25;
//        int             n_steps = 0;
//        double          delta   = 3.0;
//        Eigen::VectorXd expected(number_of_points);
//        expected << 13.6588, 12.3466, 11.0343, 9.72204, 11.3768, 13.0316, 13.0481, 13.0647, 13.0812, 13.0978, 13.1143, 13.1309, 13.1474,
//            13.164, 13.1805, 12.9996, 6.37724, 5.5777, 5.5004, 18.9982;

//        auto actual = smoothen(positions, input, f, n_steps, delta);
//        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
//    }

//    SECTION("parameter set 3") {
//        double          f       = 0.25;
//        int             n_steps = 2;
//        double          delta   = 0.0;
//        Eigen::VectorXd expected(number_of_points);
//        expected << 16.6977, 13.8516, 10.3535, 9.61856, 9.99862, 9.97798, 10.1378, 10.861, 7.20996, 7.56857, 12.1365, 16.9503, 16.7764,
//            14.0816, 13.2937, 12.9997, 7.09147, 5.93676, 5.40528, 18.998;

//        auto actual = smoothen(positions, input, f, n_steps, delta);
//        REQUIRE_THAT(actual, VectorXdIsEqualTo(expected));
//    }
//}

// TEST_CASE("Loess - Random Data") {
//    using namespace std;

//    int    points_per_period = 2000;
//    int    periods           = 10;
//    double y_scale           = 10.0;
//    int    n                 = periods * points_per_period;

//    random_device r;
//    mt19937       gen(r());

//    normal_distribution<double> d(0.0, 0.1 * y_scale);

//    auto            positions = Eigen::VectorXd::LinSpaced(n, 1000.0, 2000.0);
//    Eigen::VectorXd input(n);
//    for (int i = 0; i < n; ++i) {
//        input[i] = y_scale * sin(2.0 * M_PI * i / points_per_period) + d(gen);
//    }
//    auto output = smoothen(positions, input, 1.0 / periods / 4, 0, 10.0);

//    ofstream f("/tmp/sinus.tsv");
//    f << "x\tinput\toutput" << endl;
//    for (int i = 0; i < n; ++i) {
//        f << positions[i] << "\t" << input[i] << "\t" << output[i] << endl;
//    }
//    // TODO(sw) verify close to sinus curve
//}
