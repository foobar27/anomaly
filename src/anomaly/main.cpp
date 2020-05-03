#include "whisper/model.hpp"

#include <fstream>
#include <vector>

using namespace std;
using namespace anomaly::whisper::model;

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
        for (const auto& point : archive.m_points) {
            os << point.m_timestamp << endl;
        }
        auto time_series = archive.createTimeSeries();
        // TODO(sw) implement a range?
        for (auto point : time_series.allPoints()) {
            std::cout << "time_series[" << point.index() << "]=" << point.value() << std::endl;
        }
    }

    return 0;
}
