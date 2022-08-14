//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"

namespace vox {
/** \brief This class implements a per-triangle regular sampling for the surface
 * of 3D models.
 */
class RegularTriangleSampling {
public:
    RegularTriangleSampling();

    /** Performs the poisson sampling with the
     * respective parameters. Compare
     * http://graphics.cs.umass.edu/pubs/sa_2010.pdf
     *
     * @param numVertices number of mesh vertices
     * @param vertices vertex data of sampled data
     * @param numFaces number of faces in the mesh
     * @param faces face data of sampled mesh
     * @param maxDistance maximal distance of sampled vertices
     * @param samples vector to store the samples
     */
    static void sampleMesh(unsigned int numVertices,
                           const Vector3r *vertices,
                           unsigned int numFaces,
                           const unsigned int *faces,
                           Real maxDistance,
                           std::vector<Vector3r> &samples);

private:
    using Vector2ui = Eigen::Matrix<unsigned int, 2, 1, Eigen::DontAlign>;

    static void appendVertexSamples(unsigned int numVertices, const Vector3r *vertices, std::vector<Vector3r> &samples);
    static void appendEdgeSamples(Real d,
                                  const Vector3r *vertices,
                                  const std::vector<Vector2ui> &edges,
                                  std::vector<Vector3r> &samples,
                                  bool skipVertices = true);
    static void appendFaceSamples(Real d,
                                  const Vector3r *vertices,
                                  unsigned int numFaces,
                                  const unsigned int *faces,
                                  std::vector<Vector3r> &samples,
                                  bool skipEdges = true);

    static std::vector<Vector2ui> uniqueEdges(unsigned int numFaces, const unsigned int *faces);
};
}  // namespace vox