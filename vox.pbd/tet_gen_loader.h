//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <string>
#include <vector>

#include "vox.pbd/common.h"

namespace vox::utility {
class tet_gen_loader {
public:
    static void loadTetFile(const std::string &filename,
                            std::vector<Vector3r> &vertices,
                            std::vector<unsigned int> &tets);
    static void loadTetgenModel(const std::string &nodeFilename,
                                const std::string &eleFilename,
                                std::vector<Vector3r> &vertices,
                                std::vector<unsigned int> &tets);
    static void loadMSHModel(const std::string &mshFilename,
                             std::vector<Vector3r> &vertices,
                             std::vector<unsigned int> &tets);
};
}  // namespace vox::utility