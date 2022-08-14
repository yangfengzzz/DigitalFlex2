//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.editor/simulator_base.h"
#include "vox.sph/boundary_simulator.h"
#include "vox.sph/position_based_dynamics_wrapper/pbd_wrapper.h"

namespace vox {
class PBDBoundarySimulator : public BoundarySimulator {
protected:
    PBDWrapper *m_pbdWrapper;
    SimulatorBase *m_base;

public:
    PBDBoundarySimulator(SimulatorBase *base);
    ~PBDBoundarySimulator() override;

    void init() override;
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    void deferredInit() override;
    void timeStep() override;
    void initBoundaryData() override;
    void reset() override;

    PBDWrapper *getPBDWrapper() { return m_pbdWrapper; }
};
}  // namespace vox
