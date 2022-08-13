//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.pbd/tet_gen_loader.h"

#include <fstream>
#include <sstream>

#include "vox.base/logging.h"

using namespace vox::utility;
using namespace std;

// Call this function to load a model from a *.tet file
void tet_gen_loader::loadTetFile(const std::string &filename,
                                 std::vector<Vector3r> &vertices,
                                 std::vector<unsigned int> &tets) {
    LOGI("Loading {}", filename)

    // variables
    size_t i, num_materials, num_vertices, num_tetras, num_triangles;
    Real value;
    string line, label;
    stringstream sStream;
    // try to open the file
    ifstream fin(filename.c_str());
    if (!fin) {
        LOGE("{} file not found.", filename)
        return;
    }

    // load tet version 1.2
    getline(fin, line);
    sStream << line;
    sStream >> label;  // tet
    sStream >> label;  // version
    sStream >> value;
    sStream.clear();
    // load number of materials
    getline(fin, line);  // num_materials x
    sStream << line;
    sStream >> label;
    sStream >> num_materials;
    sStream.clear();
    // load number of vertices
    getline(fin, line);  // num_vertices x
    sStream << line;
    sStream >> label;
    sStream >> num_vertices;
    sStream.clear();
    // reverse the order of the vertices
    vertices.resize(num_vertices);
    //// read number of tetraeders
    getline(fin, line);  // num_tetras x
    sStream << line;
    sStream >> label;
    sStream >> num_tetras;
    sStream.clear();
    tets.resize(4 * num_tetras);

    // read number of triangles
    getline(fin, line);  // num_triangles x
    sStream << line;
    sStream >> label;
    sStream >> num_triangles;
    sStream.clear();

    // skip materials
    getline(fin, line);
    for (i = 0; i < num_materials; ++i) getline(fin, line);
    // skip the VERTICES label
    getline(fin, line);

    // read the vertices
    for (i = 0; i < num_vertices; ++i) {
        Real x, y, z;
        getline(fin, line);
        sStream << line;
        sStream >> x >> y >> z;
        getline(sStream, line);
        sStream.clear();

        vertices[i] = Vector3r(x, y, z);
    }

    // skip TETRAS label
    getline(fin, line);
    // read tets
    for (i = 0; i < num_tetras; ++i) {
        unsigned int tet[4];
        unsigned int m;
        getline(fin, line);
        sStream << line;
        sStream >> tet[0] >> tet[1] >> tet[2] >> tet[3] >> m;

        // for (unsigned int j = 0; j < 4; j++)
        //	tet[j]--;

        getline(sStream, line);
        sStream.clear();

        for (int j = 0; j < 4; j++) tets[4 * i + j] = tet[j];
    }
    // close file
    fin.close();

    LOGI("Number of tets: {}", num_tetras)
    LOGI("Number of vertices: {}", num_vertices)
}

void tet_gen_loader::loadTetgenModel(const std::string &nodeFilename,
                                     const std::string &eleFilename,
                                     std::vector<Vector3r> &vertices,
                                     std::vector<unsigned int> &tets) {
    LOGI("Loading {}", nodeFilename)
    LOGI("Loading {}", eleFilename)

    // variables
    size_t i, num_vertices, num_tetras;
    string nodeLine, eleLine, label;
    stringstream sStream;
    // try to open the file
    ifstream finNode(nodeFilename.c_str());
    ifstream finEle(eleFilename.c_str());
    if (!finNode) {
        LOGE("{} file not found.", nodeFilename)
        return;
    }
    if (!finEle) {
        LOGE("{} file not found.", eleFilename)
        return;
    }

    // get num vertices
    getline(finNode, nodeLine);
    sStream << nodeLine;
    sStream >> num_vertices;
    sStream >> label;  // 3
    sStream >> label;  // 0
    sStream >> label;  // 0
    sStream.clear();

    // get num tetras
    getline(finEle, eleLine);
    sStream << eleLine;
    sStream >> num_tetras;
    sStream >> label;  // 4
    sStream >> label;  // 0
    sStream >> label;  // 0
    sStream.clear();

    vertices.resize(num_vertices);
    tets.resize(4u * num_tetras);

    // read vertices
    for (i = 0; i < num_vertices; ++i) {
        unsigned nodeInd;
        Real x, y, z;
        getline(finNode, nodeLine);
        sStream << nodeLine;
        sStream >> nodeInd >> x >> y >> z;
        getline(sStream, nodeLine);
        sStream.clear();

        vertices[i] = Vector3r(x, y, z);
    }

    // read tetrahedra
    for (i = 0; i < num_tetras; ++i) {
        unsigned eleInd;
        // unsigned int tet[4];
        getline(finEle, eleLine);
        sStream << eleLine;
        sStream >> eleInd >> tets[4 * i + 0] >> tets[4 * i + 1] >> tets[4 * i + 2] >> tets[4 * i + 3];

        getline(sStream, eleLine);
        sStream.clear();
    }
    // close file
    finNode.close();
    finEle.close();

    LOGI("Number of tets: {}", num_tetras)
    LOGI("Number of vertices: {}", num_vertices)
}

void tet_gen_loader::loadMSHModel(const std::string &mshFilename,
                                  std::vector<Vector3r> &vertices,
                                  std::vector<unsigned int> &tets) {
    LOGI("Loading {}", mshFilename)

    // variables
    size_t i, num_vertices, num_tetras;
    string line, label;
    stringstream sStream;
    // try to open the file
    ifstream mshStream(mshFilename.c_str());
    if (!mshStream) {
        LOGE("{} file not found.", mshFilename)
        return;
    }

    // get num vertices
    getline(mshStream, line);
    getline(mshStream, line);
    sStream << line;
    sStream >> num_vertices;
    sStream.clear();

    vertices.resize(num_vertices);

    // read vertices
    for (i = 0; i < num_vertices; ++i) {
        unsigned nodeInd;
        Real x, y, z;
        getline(mshStream, line);
        sStream << line;
        sStream >> nodeInd >> x >> y >> z;
        getline(sStream, line);
        sStream.clear();

        vertices[i] = Vector3r(x, y, z);
    }

    // get num tetras
    getline(mshStream, line);
    getline(mshStream, line);
    getline(mshStream, line);
    sStream << line;
    sStream >> num_tetras;
    sStream.clear();

    tets.resize(4u * num_tetras);

    // read tetrahedra
    for (i = 0; i < num_tetras; ++i) {
        unsigned eleInd;
        // unsigned int tet[4];
        getline(mshStream, line);
        sStream << line;
        sStream >> eleInd >> tets[4 * i + 0] >> tets[4 * i + 1] >> tets[4 * i + 2] >> tets[4 * i + 3];

        --tets[4 * i + 0];
        --tets[4 * i + 1];
        --tets[4 * i + 2];
        --tets[4 * i + 3];

        getline(sStream, line);
        sStream.clear();
    }
    // close file
    mshStream.close();

    LOGI("Number of tets: {}", num_tetras)
    LOGI("Number of vertices: {}", num_vertices)
}
