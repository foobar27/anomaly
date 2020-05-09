#include "whisper/model.hpp"

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

constexpr double rate = 0.1;

class MedianEstimator { // TODO(sw) could be a boost accumulator
public:
    void operator()(double value) {
        if (m_first) {
            m_estimate = value;
            m_first    = false;
        } else {
            if (value > m_estimate) {
                m_estimate += rate;
            } else if (value < m_estimate) {
                m_estimate -= rate;
            }
        }
    }

    double operator*() const {
        return m_estimate;
    }

private:
    bool   m_first = true;
    double m_estimate;
};

double determineInterquantileRange(const TimeSeries& time_series) {
    using namespace boost::accumulators;
    using Accumulator = accumulator_set<double, stats<tag::extended_p_square>>;

    constexpr boost::array<double, 3> probs = {0.25, 0.75};
    Accumulator                       acc(tag::extended_p_square::probabilities = probs);
    time_series.template accumulateNonNull(acc);
    auto p25 = extended_p_square(acc)[0];
    auto p75 = extended_p_square(acc)[1];
    return p75 - p25;
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
        std::ofstream os("/tmp/percent-recent.tsv");
        auto          archive_info = archive_infos[0];
        vector<Point> points{};

        Archive archive{archive_info};
        archive.readPoints(is);
        auto time_series = archive.createTimeSeries();

        using Accumulator = boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean>>;
        int                 number_of_bins = 1440;
        vector<Accumulator> bins(number_of_bins);
        for (auto point : time_series.allPoints()) {
            auto bin = (point.m_timestamp / time_series.getConfig().m_stepInSeconds) % number_of_bins;
            if (point.m_value) {
                bins[bin](*point.m_value);
            }
        }

        auto range_in_points  = determineInterquantileRange(time_series);
        auto range_in_seconds = range_in_points * time_series.getConfig().m_stepInSeconds;
        cout << "Learning rate 1 would mean " << range_in_seconds << " seconds to traverse the inter-quantile range" << endl;

        // TODO(sw) implement a range?
        MedianEstimator median;
        MedianEstimator median_deviation;
        os << "timestamp\tvalue\tmedian\tmedian_deviation\tseason" << endl;
        for (auto point : time_series.allPoints()) {
            if (point.m_value) {
                double val    = *point.m_value;
                auto   bin    = (point.m_timestamp / time_series.getConfig().m_stepInSeconds) % number_of_bins;
                auto   season = boost::accumulators::mean(bins[bin]);

                auto noSeason = val; // TODO minus season
                median(noSeason);
                double deviation = abs(noSeason - *median);
                median_deviation(deviation);
                os << point.m_timestamp << "\t" << *point.m_value << "\t" << *median << "\t" << *median_deviation << "\t" << season << endl;
            }
        }
    }

    return 0;
}
