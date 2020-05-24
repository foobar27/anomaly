#include "timeseries.hpp"

namespace anomaly::core::timeseries {

TimeSeries::TimeSeries(const TimeSeriesConfiguration& config)
    : m_config(config)
    , m_points(config.m_numberOfPoints) {
    uint32_t ts = m_config.m_initialTimestamp;
    for (size_t i = 0; i < config.m_numberOfPoints; ++i) {
        m_points[i].m_timestamp = ts;
        ts += m_config.m_stepInSeconds;
    }
}

Eigen::VectorXd TimeSeries::toVectorXdFillNulls() const {
    // TODO(sw) what if points are missing at the beginning?
    Eigen::VectorXd result(m_points.size());
    size_t          i        = 0;
    double          previous = 0.0;
    for (const auto& point : m_points) {
        double value = previous;
        if (point.m_value)
            value = *point.m_value;
        result[i] = value;
        ++i;
    }
    return result;
}

TimeSeries derivative(const TimeSeries& time_series) {
    // TODO(sw) return a timeseries instead
    TimeSeries           output(time_series.getConfig());
    std::optional<Point> previous{};
    size_t               idx = 0;
    for (auto point : time_series.allPoints()) {
        if (point.m_value) {
            if (previous) {
                double val = (*point.m_value - *previous->m_value) / double(point.m_timestamp - previous->m_timestamp);
                // if (val < 10000) // TODO why does this happen?
                output[idx] = {point.m_timestamp, val};
            }
            previous = point;
        }
        ++idx;
    }
    return output;
}

} // end namespace anomaly::core::timeseries
