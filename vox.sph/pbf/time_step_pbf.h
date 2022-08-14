//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/pbf/simulation_data_pbf.h"
#include "vox.sph/sph_kernels.h"
#include "vox.sph/time_step.h"

namespace vox {
class SimulationDataPBF;

/** \brief This class implements the position-based fluids approach introduced
 * by Macklin and Mueller [MM13,BMO+14,BMM15].
 *
 * References:
 * - [MM13] Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1-104:12, July 2013.
 * URL: http://doi.acm.org/10.1145/2461912.2461984
 * - [BMO+14] Jan Bender, Matthias Müller, Miguel A. Otaduy, Matthias Teschner, and Miles Macklin. A survey on
 * position-based simulation methods in computer graphics. Computer Graphics Forum, 33(6):228-251, 2014. URL:
 * http://dx.doi.org/10.1111/cgf.12346
 * - [BMM15] Jan Bender, Matthias Müller, and Miles Macklin. Position-based simulation methods in computer graphics. In
 * EUROGRAPHICS 2015 Tutorials. Eurographics Association, 2015. URL: http://dx.doi.org/10.2312/egt.20151045
 */
class TimeStepPBF : public TimeStep {
protected:
    SimulationDataPBF m_simulationData;
    unsigned int m_counter;
    int m_velocityUpdateMethod;

    /** Perform a position-based correction step for the following density constraint:\n
     *  \f$C(\mathbf{x}) = \left (\frac{\rho_i}{\rho_0} - 1 \right )= 0\f$\n
     */
    void pressureSolve();
    void pressureSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();

    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

    void initParameters() override;

public:
    static int VELOCITY_UPDATE_METHOD;
    static int ENUM_PBF_FIRST_ORDER;
    static int ENUM_PBF_SECOND_ORDER;

    /** Initialize the simulation data required for this method. */
    TimeStepPBF();
    ~TimeStepPBF() override;

    /** Perform a simulation step. */
    void step() override;

    /** Reset the simulation method. */
    void reset() override;
    void resize() override;
};
}  // namespace vox
