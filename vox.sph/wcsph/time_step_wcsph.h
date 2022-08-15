//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/sph_kernels.h"
#include "vox.sph/time_step.h"
#include "vox.sph/wcsph/simulation_data_wcsph.h"

namespace vox {
class SimulationDataWCSPH;

/** \brief This class implements the Weakly Compressible SPH for Free Surface Flows approach introduced
 * by Becker and Teschner [BT07].
 *
 * References:
 * - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM
 * SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland,
 * 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
 */
class TimeStepWCSPH : public TimeStep {
protected:
    Real m_stiffness;
    Real m_exponent;

    SimulationDataWCSPH m_simulationData;
    unsigned int m_counter;

    /** Determine the pressure accelerations when the pressure is already known. */
    void computePressureAccels(unsigned int fluidModelIndex);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();

    void emittedParticles(FluidModel *model, unsigned int startIndex) override;
    void initParameters() override;

public:
    static int STIFFNESS;
    static int EXPONENT;

    TimeStepWCSPH();
    ~TimeStepWCSPH() override;

    void step() override;
    void reset() override;
    void resize() override;
};
}  // namespace vox
