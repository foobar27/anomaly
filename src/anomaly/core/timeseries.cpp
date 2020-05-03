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

} // end namespace anomaly::core::timeseries
