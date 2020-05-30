#pragma once

#include <boost/hana.hpp>
#include <fstream>

namespace anomaly::io::tsv {
template <typename RowType>
struct TsvOutputFile {
    explicit TsvOutputFile(const char* filename)
        : m_ofs(filename) {
        using namespace boost;
        RowType row{};
        bool    first = true;
        hana::for_each(row, [this, &first](auto pair) {
            if (!first) {
                m_ofs << "\t";
            }
            m_ofs << hana::to<char const*>(hana::first(pair));
            first = false;
        });
        m_ofs << std::endl;
    }

    TsvOutputFile& operator<<(const RowType& row) {
        using namespace boost;
        bool first = true;
        hana::for_each(row, [this, &first](auto pair) {
            if (!first) {
                m_ofs << "\t";
            }
            m_ofs << hana::second(pair);
            first = false;
        });
        m_ofs << std::endl;
        return *this;
    }

private:
    std::ofstream m_ofs;
};
} // end namespace anomaly::io:tsv;

#define ANOMALY_TSV_FORMAT(name, fields...) struct name { \
    BOOST_HANA_DEFINE_STRUCT(name, fields); \
};
