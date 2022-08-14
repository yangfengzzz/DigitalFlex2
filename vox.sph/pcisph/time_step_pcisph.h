//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/pcisph/simulation_data_pcisph.h"
#include "vox.sph/time_step.h"

namespace vox {
class SimulationDataPCISPH;

/** \brief This class implements the Predictive-corrective Incompressible SPH approach introduced
 * by Solenthaler and Pajarola [SP09].
 *
 * References:
 * - [SP09] B. Solenthaler and R. Pajarola. Predictive-corrective incompressible SPH. ACM Trans. Graph.,
 * 28(3):40:1-40:6, July 2009. URL: http://doi.acm.org/10.1145/1531326.1531346
 */
class TimeStepPCISPH : public TimeStep {
protected:
    SimulationDataPCISPH m_simulationData;
    unsigned int m_counter;

    void pressureSolve();
    void pressureSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();

    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

public:
    TimeStepPCISPH();
    ~TimeStepPCISPH() override;

    void step() override;
    void reset() override;
    void resize() override;
};
}  // namespace vox
