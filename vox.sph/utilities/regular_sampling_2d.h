//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "vox.sph/common.h"

namespace vox {
/** \brief This class implements a per-triangle regular sampling for the surface
 * of 3D models.
 */
class RegularSampling2D {
public:
    RegularSampling2D();

    /** Performs the poisson sampling with the
     * respective parameters. Compare
     * http://graphics.cs.umass.edu/pubs/sa_2010.pdf
     *
     * @param rotation rotation of the mesh
     * @param translation translation of the mesh
     * @param numVertices number of mesh vertices
     * @param vertices vertex data of sampled data
     * @param numFaces number of faces in the mesh
     * @param faces face data of sampled mesh
     * @param maxDistance maximal distance of sampled vertices
     * @param samples vector to store the samples
     */
    static void sampleMesh(const Matrix3r &rotation,
                           const Vector3r &translation,
                           unsigned numVertices,
                           const Vector3r *vertices,
                           unsigned int numFaces,
                           const unsigned int *faces,
                           Real maxDistance,
                           std::vector<Vector3r> &samples);
};
}  // namespace vox