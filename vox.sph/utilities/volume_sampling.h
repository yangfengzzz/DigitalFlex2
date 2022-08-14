//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"

namespace vox::utility {
/** \brief This class implements a volume sampling of 3D models.
 */
class VolumeSampling {
public:
    /** Performs the volume sampling with the
     * respective parameters.
     *
     * @param numVertices number of vertices
     * @param vertices vertex data
     * @param numFaces number of faces
     * @param faces index list of faces
     * @param radius radius of sampled particles
     * @param region defines a subregion of the mesh to be sampled (nullptr if not used)
     * @param resolution resolution of the used SDF
     * @param invert defines if the mesh should be inverted and the outside is sampled
     * @param sampleMode 0=regular, 1=almost dense, 2=dense
     * @param samples sampled vertices that will be returned
     */
    static void sampleMesh(unsigned int numVertices,
                           const Vector3r *vertices,
                           unsigned int numFaces,
                           const unsigned int *faces,
                           Real radius,
                           const AlignedBox3r *region,
                           const std::array<unsigned int, 3> &resolution,
                           bool invert,
                           unsigned int sampleMode,
                           std::vector<Vector3r> &samples);
};
}  // namespace vox::utility