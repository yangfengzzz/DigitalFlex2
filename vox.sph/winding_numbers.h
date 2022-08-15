//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/triangle_mesh.h"

namespace vox::utility {
class WindingNumbers {
public:
    /** Determine the winding number for a point p and a triangle abc.
     */
    static Real computeGeneralizedWindingNumber(const Vector3r& p,
                                                const Vector3r& a,
                                                const Vector3r& b,
                                                const Vector3r& c);

    /** Determine the winding number of a point p in a triangle mesh.
     */
    static Real computeGeneralizedWindingNumber(const Vector3r& p, const TriangleMesh& mesh);
};
}  // namespace vox::utility