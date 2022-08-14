//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/iisph/simulation_data_iisph.h"
#include "vox.sph/sph_kernels.h"
#include "vox.sph/time_step.h"

namespace vox {
class SimulationDataIISPH;

/** \brief This class implements the Implicit Incompressible SPH approach introduced
 * by Ihmsen et al. [ICS+14].
 *
 * References:
 * - [ICS+14] Markus Ihmsen, Jens Cornelis, Barbara Solenthaler, Christopher Horvath, and Matthias Teschner. Implicit
 * incompressible SPH. IEEE Transactions on Visualization and Computer Graphics, 20(3):426-435, March 2014. URL:
 * http://dx.doi.org/10.1109/TVCG.2013.105
 */
class TimeStepIISPH : public TimeStep {
protected:
    SimulationDataIISPH m_simulationData;
    unsigned int m_counter;

    void predictAdvection(unsigned int fluidModelIndex);
    void pressureSolve();
    void pressureSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);
    void integration(unsigned int fluidModelIndex);

    /** Determine the pressure accelerations when the pressure is already known. */
    void computePressureAccels(unsigned int fluidModelIndex);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();

    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

public:
    TimeStepIISPH();
    ~TimeStepIISPH() override;

    void step() override;
    void reset() override;
    void resize() override;

    const SimulationDataIISPH &getSimulationData() { return m_simulationData; };
};
}  // namespace vox