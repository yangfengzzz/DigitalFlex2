//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/parameter_object.h"
#include "vox.pbd/common.h"
#include "vox.pbd/simulation_model.h"
#include "vox.pbd/time_step.h"

namespace vox {
/** \brief Class to manage the current simulation time and the time step size.
 * This class is a singleton.
 */
class simulation : public ParameterObject {
public:
    static int GRAVITATION;

protected:
    SimulationModel *m_model{};
    TimeStep *m_timeStep;
    Vector3r m_gravitation;

    void initParameters() override;

private:
    static simulation *current;

public:
    simulation();
    simulation(const simulation &) = delete;
    simulation &operator=(const simulation &) = delete;
    ~simulation() override;

    void init();
    void reset();

    // Singleton
    static simulation *getCurrent();
    static void setCurrent(simulation *tm);
    static bool hasCurrent();

    SimulationModel *getModel() { return m_model; }
    void setModel(SimulationModel *model) { m_model = model; }

    TimeStep *getTimeStep() { return m_timeStep; }
    void setTimeStep(TimeStep *ts) { m_timeStep = ts; }
};
}  // namespace vox