#pragma once

#include <optional>
#include <vector>
#include <boost/accumulators/accumulators.hpp>

namespace anomaly::core::timeseries {

struct TimeSeriesConfiguration {
    size_t   m_numberOfPoints;
    uint32_t m_initialTimestamp;
    uint32_t m_finalTimestamp;
    size_t   m_stepInSeconds;
};

struct Point {
    uint32_t              m_timestamp;
    std::optional<double> m_value;
};

// TODO(sw) rotating
class TimeSeries {

public:
    TimeSeries(const TimeSeriesConfiguration& config);

    const TimeSeriesConfiguration& getConfig() const {
        return m_config;
    }

    Point& operator[](size_t i) {
        return m_points[i];
    }

    const Point& operator[](size_t i) const {
        return m_points[i];
    }

    auto allPoints() {
        return m_points;
    }

    auto allPoints() const {
        return m_points;
    }

    // TODO inline this once we have proper ranges
    template <typename Accumulator>
    void accumulateNonNull(Accumulator& acc) const {
        // TODO simplify with ranges and std::for_each(..., acc)
        for (auto& point : m_points) {
            if (point.m_value) {
                acc(*point.m_value);
            }
        }
    }

private:
    TimeSeriesConfiguration m_config;
    std::vector<Point>      m_points;
};

} // end namespace anomaly::core::timeseries
