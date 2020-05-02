#include "model.hpp"

#include <boost/endian.hpp>
#include <boost/integer.hpp>

namespace anomaly::whisper::model {

std::ostream & operator << (std::ostream &out, const DebugWrapper<AggregationType> & t) {
    switch (t.m_value) {
    case AggregationType::average:
        out << "average";
        break;
    case AggregationType::sum:
        out << "sum";
        break;
    case AggregationType::last:
        out << "last";
        break;
    case AggregationType::max:
        out << "max";
        break;
    case AggregationType::min:
        out << "min";
        break;
    }
    return out;
}

//
// Read/Write big endian values
//

//template<typename T>
//static void writeBig(std::ostream & out, const T &value) {
//    auto raw_value = boost::endian::native_to_big(value);
//    out.write(reinterpret_cast<const char *>(&raw_value), sizeof(raw_value));
//}

// Workaround for flooat and double, which boost refuses to convert.
template<typename T>
using raw_type = typename boost::uint_t<8 * sizeof(T)>::exact;

template<typename T>
static void writeBig(std::ostream & out, const T &value) {
    raw_type<T> raw_value {};
    std::memcpy(&raw_value, &value, sizeof(value));
    boost::endian::native_to_big_inplace(raw_value);
    // TODO(sw) change to bit_cast
    out.write(reinterpret_cast<const char *>(&raw_value), sizeof(raw_value));
}

//template<typename T>
//static void readBig(std::istream & in, T &value) {
//    in.read(reinterpret_cast<char *>(&value), sizeof(value));
//    boost::endian::big_to_native_inplace(value);
//}

// TODO(sw) somehow get rid of this copy step
template<typename T>
static void readBig(std::istream & in, T &value) {
    raw_type<T> raw_value {};
    in.read(reinterpret_cast<char *>(&raw_value), sizeof(raw_value));
    boost::endian::big_to_native_inplace(raw_value);
    std::memcpy(&value, &raw_value, sizeof(value));
}

//
// MetaData
//

std::ostream & operator << (std::ostream &out, const MetaData &m) {
    writeBig(out, m.m_aggregationType);
    writeBig(out, m.m_maxRetention);
    writeBig(out, m.m_xFilesFactor);
    writeBig(out, m.m_archiveCount);
    return out;
}

std::istream & operator >> (std::istream &in, MetaData & m) {
    readBig(in, m.m_aggregationType);
    readBig(in, m.m_maxRetention);
    readBig(in, m.m_xFilesFactor);
    readBig(in, m.m_archiveCount);
    return in;
}

std::ostream & operator << (std::ostream &out, const DebugWrapper<MetaData> & m) {
    out << "MetaData{aggregationType: " << debug(m.m_value.m_aggregationType)
        << ", maxRetention: " << m.m_value.m_maxRetention
        << ", xFilesFactor: " << m.m_value.m_xFilesFactor
        << ", archiveCount: " << m.m_value.m_archiveCount
        << "}";
    return out;
}

//
// ArchiveInfo
//

std::ostream & operator << (std::ostream &out, const ArchiveInfo &a) {
    writeBig(out, a.m_offset);
    writeBig(out, a.m_secondsPerPoint);
    writeBig(out, a.m_numberOfPoints);
    return out;
}

std::istream & operator >> (std::istream &in, ArchiveInfo &a) {
    readBig(in, a.m_offset);
    readBig(in, a.m_secondsPerPoint);
    readBig(in, a.m_numberOfPoints);
    return in;
}

std::ostream & operator << (std::ostream &out, const DebugWrapper<ArchiveInfo> &a) {
    out << "ArchiveInfo{offset: " << a.m_value.m_offset
        << ", secondsPerPoint: " << a.m_value.m_secondsPerPoint
        << ", numberOfPoints: " << a.m_value.m_numberOfPoints
        << "}";
    return out;
}

//
// Point
//

std::ostream & operator << (std::ostream &out, const Point &p) {
    writeBig(out, p.m_timestamp);
    writeBig(out, p.m_value);
    return out;
}

std::istream & operator >> (std::istream &in, Point &p) {
    readBig(in, p.m_timestamp);
    readBig(in, p.m_value);
    return in;
}

std::ostream & operator << (std::ostream &out, const DebugWrapper<Point> &p) {
    out << "Point{timestamp: " << p.m_value.m_timestamp
        << ", value: " << p.m_value.m_value
        << "}";
    return out;
}

}
