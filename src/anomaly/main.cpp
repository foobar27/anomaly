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

    // TODO(sw) how to iterate in order: rotating files, reverse archives?
    for (auto & archive_info : archive_infos) {
        for (uint i = 0; i < archive_info.m_numberOfPoints; ++i) {
            Point point {};
            is >> point;
            //std::cout << "point[" << i << "] = ";
            //cout << debug(point) << endl;
        }
    }

    return 0;
}
