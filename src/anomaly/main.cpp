#include "whisper/model.hpp"

#include <fstream>
#include <vector>

using namespace std;
using namespace anomaly::whisper::model;

int main() {
    ifstream is ("/home/sebastien/percent.wsp", ifstream::binary);
    MetaData meta_data {};
    is >> meta_data;

    cout << debug(meta_data) << endl;

    vector<ArchiveInfo> archive_infos {};
    for (uint i = 0; i < meta_data.m_archiveCount; ++i) {
        ArchiveInfo archive_info {};
        is >> archive_info;
        archive_infos.push_back(archive_info);
        cout << debug(archive_info)
             << " retention=" << archive_info.retention()
             << " size=" << archive_info.size()
             <<  endl;
    }

    //for (auto & archive_info : archive_infos) {
    {
        std::ofstream os("/tmp/percent-recent.tsv");
        auto archive_info = archive_infos[0];
        vector<Point> points {};

        Archive archive { archive_info };
        archive.readPoints(is);
        uint32_t ts {};
        for (const auto & point : archive.m_points) {
            auto delta = (point.m_timestamp - ts);
            if (delta != 0) {
                std::cout << delta << endl;
            }
            ts = point.m_timestamp;
        }
    }

    return 0;
}
