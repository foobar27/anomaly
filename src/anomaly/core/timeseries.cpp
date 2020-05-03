#include "timeseries.hpp"

namespace anomaly::core::timeseries {

TimeSeries::TimeSeries(const TimeSeriesConfiguration& config)
    : m_config(config)
    , m_values(config.m_numberOfPoints)
    , m_initialized(config.m_numberOfPoints, 0) { }

bool TimeSeries::isInitialized(size_t idx) const {
    return m_initialized[idx];
}

std::optional<double> TimeSeries::getAtIndex(size_t idx) const {
    if (m_initialized[idx]) {
        return m_values[idx];
    }
    return m_values[idx];
}

void TimeSeries::setAtIndex(size_t idx, double value) {
    m_initialized[idx] = true;
    m_values[idx]      = value;
}

void TimeSeries::resetAtIndex(size_t idx) {
    m_initialized[idx] = false;
}

} // end namespace anomaly::core::timeseries
