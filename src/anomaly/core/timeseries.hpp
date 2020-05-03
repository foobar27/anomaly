#pragma once

#include <optional>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/indexed.hpp>

namespace anomaly::core::timeseries {

struct TimeSeriesConfiguration {
    size_t   m_numberOfPoints;
    uint32_t m_initialTimestamp;
    uint32_t m_finalTimestamp;
    size_t   m_stepInSeconds;
};

using ValueVector        = std::vector<double>;
using ValueConstIterator = typename ValueVector::const_iterator;

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

    auto allPoints() const {
        using namespace boost::adaptors;
        return m_values | indexed(0);
    }

private:
    TimeSeriesConfiguration m_config;
    ValueVector             m_values;
    boost::dynamic_bitset<> m_initialized;
};

} // end namespace anomaly::core::timeseries
