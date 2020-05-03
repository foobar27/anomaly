#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

#include "../core/timeseries.hpp"

namespace anomaly::whisper::model {

template <typename T>
struct DebugWrapper {
    const T& m_value;
};

template <typename T>
DebugWrapper<T> debug(const T& value) {
    return {value};
}

enum class AggregationType : uint32_t { average = 1, sum = 2, last = 3, max = 4, min = 5 };

std::ostream& operator<<(std::ostream& out, const DebugWrapper<AggregationType>&);

struct MetaData {
    AggregationType m_aggregationType;
    uint32_t        m_maxRetention;
    float           m_xFilesFactor;
    uint32_t        m_archiveCount;
};

std::ostream& operator<<(std::ostream&, const MetaData&);
std::istream& operator>>(std::istream&, MetaData&);
std::ostream& operator<<(std::ostream& out, const DebugWrapper<MetaData>& m);

struct ArchiveInfo {
    uint32_t m_offset;
    uint32_t m_secondsPerPoint;
    uint32_t m_numberOfPoints;

    uint32_t retention() const {
        return m_secondsPerPoint * m_numberOfPoints;
    }

    uint32_t size() const {
        constexpr auto pointSize = 4 + 8; // uint32_t(4) + double(8)
        return m_numberOfPoints * pointSize;
    }
};

std::ostream& operator<<(std::ostream&, const anomaly::whisper::model::ArchiveInfo&);
std::istream& operator>>(std::istream&, anomaly::whisper::model::ArchiveInfo&);
std::ostream& operator<<(std::ostream& out, const anomaly::whisper::model::DebugWrapper<anomaly::whisper::model::ArchiveInfo>& a);

struct Point {
    uint32_t m_timestamp;
    double   m_value;
};

std::ostream& operator<<(std::ostream&, const Point&);
std::istream& operator>>(std::istream&, Point&);
std::ostream& operator<<(std::ostream& out, const DebugWrapper<Point>& p);

struct Archive {
    ArchiveInfo        m_info;
    std::vector<Point> m_points;

    Archive(ArchiveInfo& info)
        : m_info{info}
        , m_points{m_info.m_numberOfPoints} { }

    void readPoints(std::istream&);

    core::timeseries::TimeSeries createTimeSeries() const;
};

} // end namespace anomaly::whisper::model
