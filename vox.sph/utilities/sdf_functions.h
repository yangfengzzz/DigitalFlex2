//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/discrete_grid/cubic_lagrange_discrete_grid.h"

namespace vox::utility {
/** \brief Functions for generating and querying an SDF.
 */
class SDFFunctions {
public:
    /** Generate SDF from mesh.
     */
    static CubicLagrangeDiscreteGrid *generateSDF(unsigned int numVertices,
                                                  const Vector3r *vertices,
                                                  unsigned int numFaces,
                                                  const unsigned int *faces,
                                                  const AlignedBox3r &bbox,
                                                  const std::array<unsigned int, 3> &resolution,
                                                  bool invert = false);

    /** Compute the bounding box of a mesh.
     */
    static AlignedBox3r computeBoundingBox(unsigned int numVertices, const Vector3r *vertices);

    /** Determine distance of a point x to the surface represented by the SDF and corresponding surface normal and
     * next point on the surface.
     */
    static double distance(CubicLagrangeDiscreteGrid *sdf,
                           const Vector3r &x,
                           Real thickness,
                           Vector3r &normal,
                           Vector3r &nextSurfacePoint);

    /** Determine distance of a point x to the surface represented by the SDF.
     */
    static double distance(CubicLagrangeDiscreteGrid *sdf, const Vector3r &x, Real thickness);
};
}  // namespace vox::utility