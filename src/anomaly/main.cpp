#include "whisper/model.hpp"
#include "anomaly/core/stl.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/array.hpp>

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace anomaly::whisper::model;
using anomaly::core::timeseries::TimeSeries;

class MedianEstimator { // TODO(sw) could be a boost accumulator
public:
    explicit MedianEstimator(double learningRate)
        : m_learningRate{learningRate}
        , m_first{true}
        , m_estimate{0} { }

    void operator()(double value) {
        if (m_first) {
            m_estimate = value;
            m_first    = false;
        } else {
            if (value > m_estimate) {
                m_estimate += m_learningRate;
            } else if (value < m_estimate) {
                m_estimate -= m_learningRate;
            }
        }
    }

    double operator*() const {
        return m_estimate;
    }

private:
    const double m_learningRate;
    bool         m_first;
    double       m_estimate;
};

static double determineInterquantileRange(const TimeSeries& time_series) {
    using namespace boost::accumulators;
    using Accumulator = accumulator_set<double, stats<tag::extended_p_square>>;

    constexpr boost::array<double, 3> probs = {0.25, 0.75};
    Accumulator                       acc(tag::extended_p_square::probabilities = probs);
    time_series.template accumulateNonNull(acc);
    auto p25 = extended_p_square(acc)[0];
    auto p75 = extended_p_square(acc)[1];
    return p75 - p25;
}

// TODO(sw) this should be a cyclic TimeSeries (less storage)
static TimeSeries determineSeasonalComponent(const TimeSeries& input) {
    using Accumulator = boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean>>;
    auto                step_in_seconds = input.getConfig().m_stepInSeconds;
    size_t              number_of_bins  = 1440; // TODO(sw) deduce
    vector<Accumulator> bins(number_of_bins);
    for (auto& point : input.allPoints()) {
        auto bin = (point.m_timestamp / step_in_seconds) % number_of_bins;
        if (point.m_value) {
            bins[bin](*point.m_value);
        }
    }
    TimeSeries output{input.getConfig()};
    for (auto& point : input.allPoints()) {
        auto bin = (point.m_timestamp / step_in_seconds) % number_of_bins;
        // TODO(sw) only set if bin is not empty!
        point.m_value = boost::accumulators::mean(bins[bin]);
    }
    return output;
}

int main0() {
    using namespace anomaly::core::stl;
    ifstream is("/home/sebastien/percent.wsp", ifstream::binary);
    MetaData meta_data{};
    is >> meta_data;

    cout << debug(meta_data) << endl;

    vector<ArchiveInfo> archive_infos{};
    for (uint i = 0; i < meta_data.m_archiveCount; ++i) {
        ArchiveInfo archive_info{};
        is >> archive_info;
        archive_infos.push_back(archive_info);
        cout << debug(archive_info);
        cout << " retention=" << archive_info.retention();
        cout << " size=" << archive_info.size();
        cout << endl;
    }

    Archive archive{archive_infos[0]};
    archive.readPoints(is);
    auto time_series = archive.createTimeSeries();

    Eigen::VectorXd input = time_series.toVectorXdFillNulls();

    int  points_per_period = 1440;
    auto low_pass_length   = points_per_period;
    if (low_pass_length % 2 == 0)
        low_pass_length++;
    // smoothers: length, degree, ratio, number of steps, delta
    StlConfiguration config{points_per_period,
                            // seasonal smoother
                            {7, 1, 0.0, 0, 1},
                            // trend smoother
                            {points_per_period * 2, 1, 0.0, 0, points_per_period / 20.0},
                            // lowpass smoother
                            {low_pass_length, 1, 0.0, 0, points_per_period / 20.0},
                            // number of iterations inner loop
                            1,
                            // number of iterations outer loop
                            1};

    StlAlgorithm decomposition(config, input.size());
    {
        // Performance testing
        auto start        = std::chrono::high_resolution_clock::now();
        auto n_iterations = 10;
        for (int i = 0; i < n_iterations; ++i) {
            decomposition.stl(input);
        }
        auto finish       = std::chrono::high_resolution_clock::now();
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
        std::cout << (microseconds.count() / n_iterations) << "µs\n";
    }

    decomposition.stl(input);

    ofstream f("/tmp/percent-decomposed.tsv");
    f << "x\tinput\tseason\ttrend" << endl;
    for (int i = 0; i < input.size(); ++i) {
        f << i << "\t" << input[i] << "\t" << decomposition.season()[i] << "\t" << decomposition.trend()[i] << endl;
    }

    return 0;
}

int main2() {
    using namespace Eigen;
    using namespace std;
    VectorXd big(5);
    big << 10, 11, 12, 13, 14;
    auto segment = big.segment(1, 2);
    segment[0]   = 100;
    cout << segment.sum() << endl;
    cout << big[1] << endl;
    return 0;
}

int main() {
    ifstream is("/home/sebastien/end-offset.wsp", ifstream::binary);
    MetaData meta_data{};
    is >> meta_data;

    cout << debug(meta_data) << endl;

    vector<ArchiveInfo> archive_infos{};
    for (uint i = 0; i < meta_data.m_archiveCount; ++i) {
        ArchiveInfo archive_info{};
        is >> archive_info;
        archive_infos.push_back(archive_info);
        cout << debug(archive_info);
        cout << " retention=" << archive_info.retention();
        cout << " size=" << archive_info.size();
        cout << endl;
    }

    // for (auto & archive_info : archive_infos) {
    {
        ofstream      os("/tmp/end-offset-recent.tsv");
        auto          archive_info = archive_infos[0];
        vector<Point> points{};

        Archive archive{archive_info};
        archive.readPoints(is);
        auto time_series = archive.createTimeSeries();

        auto range_in_points  = determineInterquantileRange(time_series);
        auto range_in_seconds = range_in_points * double(time_series.getConfig().m_stepInSeconds);
        cout << "Learning rate 1 would mean " << range_in_seconds << " seconds to traverse the inter-quantile range" << endl;

        auto diff = derivative(time_series);

        auto seasonal = determineSeasonalComponent(time_series);
        // TODO(sw) implement a range?
        MedianEstimator median(0.001);
        MedianEstimator median_deviation(0.001);
        os << "timestamp\tvalue\tmedian\tmedian_deviation\tseason" << endl;
        for (auto point : diff.allPoints()) {
            double val    = *point.m_value;
            auto   season = 0.0; // seasonal_it->m_value;

            auto no_season = val;
            median(no_season);
            double deviation = abs(no_season - *median);
            median_deviation(deviation);
            os << point.m_timestamp << "\t" << val << "\t" << *median << "\t" << *median_deviation << "\t" << season << endl;
        }
    }

    return 0;
}
