//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/sph_kernels.h"
#include "vox.sph/icsph/simulation_data_icsph.h"
#include "vox.sph/time_step.h"

namespace vox {
class SimulationDataICSPH;

/** \brief This class implements the Implicit Compressible SPH approach introduced
 * by Gissler et al. [GHB+20].
 *
 * References:
 * - [GHB+20] Christoph Gissler, Andreas Henne, Stefan Band, Andreas Peer and Matthias Teschner. An Implicit
 * Compressible SPH Solver for Snow Simulation. ACM Transactions on Graphics, 39(4). URL:
 * https://doi.org/10.1145/3386569.3392431
 */
class TimeStepICSPH : public TimeStep {
protected:
    SimulationDataICSPH m_simulationData;
    Real m_lambda;
    bool m_clamping;
    const Real m_psi = 1.5;
    unsigned int m_counter;

    void computeDensityAdv(unsigned int fluidModelIndex);
    void compute_aii(unsigned int fluidModelIndex);
    void pressureSolve();
    void pressureSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);
    void integration(unsigned int fluidModelIndex);

    /** Determine the pressure accelerations when the pressure is already known. */
    void computePressureAccels(unsigned int fluidModelIndex);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();

    void initParameters() override;
    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

public:
    static int LAMBDA;
    static int PRESSURE_CLAMPING;

    TimeStepICSPH();
    ~TimeStepICSPH() override;

    void step() override;
    void reset() override;
    void resize() override;

    const SimulationDataICSPH &getSimulationData() { return m_simulationData; };
};
}  // namespace vox
