#include "whisper/model.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <boost/array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/stats.hpp>

using namespace std;
using namespace anomaly::whisper::model;
using anomaly::core::timeseries::TimeSeries;

class MedianEstimator { // TODO(sw) could be a boost accumulator
public:
    MedianEstimator(double learningRate)
        : m_learningRate(learningRate) { }

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
    bool         m_first = true;
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

// TODO this should be a cyclic TimeSeries (less storage)
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
        // TODO only set if bin is not empty!
        point.m_value = boost::accumulators::mean(bins[bin]);
    }
    return output;
}

#include <eigen3/Eigen/Dense>

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

    // for (auto & archive_info : archive_infos) {
    {
        ofstream      os("/tmp/percent-recent.tsv");
        auto          archive_info = archive_infos[0];
        vector<Point> points{};

        Archive archive{archive_info};
        archive.readPoints(is);
        auto time_series = archive.createTimeSeries();

        auto range_in_points  = determineInterquantileRange(time_series);
        auto range_in_seconds = range_in_points * double(time_series.getConfig().m_stepInSeconds);
        cout << "Learning rate 1 would mean " << range_in_seconds << " seconds to traverse the inter-quantile range" << endl;

        auto seasonal = determineSeasonalComponent(time_series);
        // TODO(sw) implement a range?
        MedianEstimator median(1.0);
        MedianEstimator median_deviation(1.0);
        os << "timestamp\tvalue\tmedian\tmedian_deviation\tseason" << endl;
        auto seasonal_it = seasonal.allPoints().begin(); // TODO(sw) iterate over two ranges
        for (auto point : time_series.allPoints()) {
            if (point.m_value) {
                double val    = *point.m_value;
                auto   season = seasonal_it->m_value;

                auto noSeason = val; // TODO(sw) minus season
                median(noSeason);
                double deviation = abs(noSeason - *median);
                median_deviation(deviation);
                os << point.m_timestamp << "\t" << *point.m_value << "\t" << *median << "\t" << *median_deviation << "\t" << *season << endl;
            }
            seasonal_it++;
        }
    }

    return 0;
}
