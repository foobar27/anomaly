#pragma once

#include <optional>
#include <vector>
#include <boost/dynamic_bitset.hpp>

namespace anomaly::core::timeseries {

struct TimeSeriesConfiguration {
    size_t   m_numberOfPoints;
    uint32_t m_initialTimestamp;
    uint32_t m_finalTimestamp;
    size_t   m_stepInSeconds;
};

using ValueVector        = std::vector<double>;
using ValueConstIterator = typename ValueVector::const_iterator;

class RawValueRange {
public:
    RawValueRange(const ValueVector& values)
        : m_values(values) { }

    ValueConstIterator begin() {
        return m_values.begin();
    }

    ValueConstIterator end() {
        return m_values.end();
    }

private:
    const ValueVector& m_values;
};

// TODO(sw) rotating
class TimeSeries {

public:
    TimeSeries(const TimeSeriesConfiguration& config);

    const TimeSeriesConfiguration& getConfig() const {
        return m_config;
    }

    bool isInitialized(size_t idx) const;

    std::optional<double> getAtIndex(size_t idx) const;

    void setAtIndex(size_t idx, double value);

    void resetAtIndex(size_t idx);

    RawValueRange rawValueRange() const {
        return {m_values};
    }

private:
    TimeSeriesConfiguration m_config;
    ValueVector             m_values;
    boost::dynamic_bitset<> m_initialized;
};

} // end namespace anomaly::core::timeseries
