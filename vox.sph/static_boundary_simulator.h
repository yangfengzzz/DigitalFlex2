//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/boundary_simulator.h"

namespace vox {
class SimulatorBase;
class TriangleMesh;

class StaticBoundarySimulator : public BoundarySimulator {
protected:
    SimulatorBase *m_base;

    static void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale);

public:
    explicit StaticBoundarySimulator(SimulatorBase *base);
    ~StaticBoundarySimulator() override;

    void initBoundaryData() override;
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    void deferredInit() override;

    void timeStep() override;
    void reset() override;
};
}  // namespace vox