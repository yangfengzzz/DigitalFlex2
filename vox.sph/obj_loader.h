//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <array>
#include <string>

#include "vox.base/logging.h"
#include "vox.base/string_tools.h"

namespace vox::utility {
/** \brief Struct to store the position and normal indices
 */
struct MeshFaceIndices {
    int posIndices[3];
    int texIndices[3];
    int normalIndices[3];
};

/** \brief Read for OBJ files.
 */
class OBJLoader {
public:
    using Vec3f = std::array<float, 3>;
    using Vec2f = std::array<float, 2>;

    /** This function loads an OBJ file.
     * Only triangulated meshes are supported.
     */

    static void loadObj(const std::string &filename,
                        std::vector<Vec3f> *x,
                        std::vector<MeshFaceIndices> *faces,
                        std::vector<Vec3f> *normals,
                        std::vector<Vec2f> *texcoords,
                        const Vec3f &scale) {
        LOGI("Loading {}", filename)

        std::ifstream filestream;
        filestream.open(filename.c_str());
        if (filestream.fail()) {
            LOGE("Failed to open file: {}", filename)
            return;
        }

        std::string line_stream;
        bool vt = false;
        bool vn = false;

        std::vector<std::string> pos_buffer;
        std::vector<std::string> f_buffer;

        while (getline(filestream, line_stream)) {
            std::stringstream str_stream(line_stream);
            std::string type_str;
            str_stream >> type_str;

            if (type_str == "v") {
                Vec3f pos;
                pos_buffer.clear();
                std::string parse_str = line_stream.substr(line_stream.find('v') + 1);
                StringTools::tokenize(parse_str, pos_buffer);
                for (unsigned int i = 0; i < 3; i++) pos[i] = stof(pos_buffer[i]) * scale[i];

                x->push_back(pos);
            } else if (type_str == "vt") {
                if (texcoords != nullptr) {
                    Vec2f tex;
                    pos_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find("vt") + 2);
                    StringTools::tokenize(parse_str, pos_buffer);
                    for (unsigned int i = 0; i < 2; i++) tex[i] = stof(pos_buffer[i]);

                    texcoords->push_back(tex);
                    vt = true;
                }
            } else if (type_str == "vn") {
                if (normals != nullptr) {
                    Vec3f nor;
                    pos_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find("vn") + 2);
                    StringTools::tokenize(parse_str, pos_buffer);
                    for (unsigned int i = 0; i < 3; i++) nor[i] = stof(pos_buffer[i]);

                    normals->push_back(nor);
                    vn = true;
                }
            } else if (type_str == "f") {
                MeshFaceIndices faceIndex{};
                if (vn && vt) {
                    f_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find('f') + 1);
                    StringTools::tokenize(parse_str, f_buffer);
                    for (int i = 0; i < 3; ++i) {
                        pos_buffer.clear();
                        StringTools::tokenize(f_buffer[i], pos_buffer, "/");
                        faceIndex.posIndices[i] = stoi(pos_buffer[0]);
                        faceIndex.texIndices[i] = stoi(pos_buffer[1]);
                        faceIndex.normalIndices[i] = stoi(pos_buffer[2]);
                    }
                } else if (vn) {
                    f_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find('f') + 1);
                    StringTools::tokenize(parse_str, f_buffer);
                    for (int i = 0; i < 3; ++i) {
                        pos_buffer.clear();
                        StringTools::tokenize(f_buffer[i], pos_buffer, "/");
                        faceIndex.posIndices[i] = stoi(pos_buffer[0]);
                        faceIndex.normalIndices[i] = stoi(pos_buffer[1]);
                    }
                } else if (vt) {
                    f_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find('f') + 1);
                    StringTools::tokenize(parse_str, f_buffer);
                    for (int i = 0; i < 3; ++i) {
                        pos_buffer.clear();
                        StringTools::tokenize(f_buffer[i], pos_buffer, "/");
                        faceIndex.posIndices[i] = stoi(pos_buffer[0]);
                        faceIndex.texIndices[i] = stoi(pos_buffer[1]);
                    }
                } else {
                    f_buffer.clear();
                    std::string parse_str = line_stream.substr(line_stream.find('f') + 1);
                    StringTools::tokenize(parse_str, f_buffer);
                    for (int i = 0; i < 3; ++i) {
                        faceIndex.posIndices[i] = stoi(f_buffer[i]);
                    }
                }
                faces->push_back(faceIndex);
            }
        }
        filestream.close();
    }
};
}  // namespace vox::utility
