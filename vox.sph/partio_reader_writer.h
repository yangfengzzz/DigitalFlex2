//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"

namespace vox::utility {
/** \brief Class for reading and writing partio files.
 */
class PartioReaderWriter {
public:
    static bool readParticles(const std::string &fileName,
                              const Vector3r &translation,
                              const Matrix3r &rotation,
                              Real scale,
                              std::vector<Vector3r> &pos,
                              std::vector<Vector3r> &vel);

    static bool readParticles(const std::string &fileName,
                              const Vector3r &translation,
                              const Matrix3r &rotation,
                              Real scale,
                              std::vector<Vector3r> &positions,
                              std::vector<Vector3r> &velocities,
                              Real &particleRadius);

    static bool readParticles(const std::string &fileName,
                              const Vector3r &translation,
                              const Matrix3r &rotation,
                              Real scale,
                              std::vector<Vector3r> &pos);

    static void writeParticles(const std::string &fileName,
                               unsigned int numParticles,
                               const Vector3r *particlePositions,
                               const Vector3r *particleVelocities,
                               Real particleRadius);
};

}  // namespace vox::utility