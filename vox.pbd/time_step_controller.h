//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.pbd/collision_detection.h"
#include "vox.pbd/simulation_model.h"
#include "vox.pbd/time_step.h"

namespace vox {
class TimeStepController : public TimeStep {
public:
    // 		static int SOLVER_ITERATIONS;
    // 		static int SOLVER_ITERATIONS_V;
    static int NUM_SUB_STEPS;
    static int MAX_ITERATIONS;
    static int MAX_ITERATIONS_V;
    static int VELOCITY_UPDATE_METHOD;

    static int ENUM_VUPDATE_FIRST_ORDER;
    static int ENUM_VUPDATE_SECOND_ORDER;

protected:
    int m_velocityUpdateMethod;
    unsigned int m_iterations;
    unsigned int m_iterationsV;
    unsigned int m_subSteps;
    unsigned int m_maxIterations;
    unsigned int m_maxIterationsV;

    void initParameters() override;

    void positionConstraintProjection(SimulationModel &model);
    void velocityConstraintProjection(SimulationModel &model);

public:
    TimeStepController();
    ~TimeStepController() override;

    void step(SimulationModel &model) override;
    void reset() override;
};
}  // namespace vox