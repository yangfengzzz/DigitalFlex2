//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"

namespace vox {
class BoundarySimulator {
public:
    BoundarySimulator() = default;
    virtual ~BoundarySimulator() = default;
    virtual void init() {}
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    virtual void deferredInit() {}
    virtual void timeStep() {}
    virtual void initBoundaryData() {}
    virtual void reset() {}

    static void updateBoundaryForces();
};
}  // namespace vox