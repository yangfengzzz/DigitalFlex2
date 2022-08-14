//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/dfsph/simulation_data_dfsph.h"
#include "vox.sph/sph_kernels.h"
#include "vox.sph/time_step.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace vox {
class SimulationDataDFSPH;

/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
 * by Bender and Koschier [BK15,BK17,KBST19].
 *
 * References:
 * - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics
 * Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL:
 * http://doi.acm.org/10.1145/2786784.2786796
 * - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on
 * Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
 * - [KBST19] Dan Koschier, Jan Bender, Barbara Solenthaler, and Matthias Teschner. Smoothed particle hydrodynamics for
 * physically-based simulation of fluids and solids. In Eurographics 2019 - Tutorials. Eurographics Association, 2019.
 * URL: https://interactivecomputergraphics.github.io/SPH-Tutorial
 */
class TimeStepDFSPH : public TimeStep {
protected:
    SimulationDataDFSPH m_simulationData;
    unsigned int m_counter;
    const Real m_eps = static_cast<Real>(1.0e-5);
    bool m_enableDivergenceSolver;
    unsigned int m_iterationsV;
    Real m_maxErrorV;
    unsigned int m_maxIterationsV;

    void computeDFSPHFactor(unsigned int fluidModelIndex);
    void pressureSolve();
    void pressureSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);
    void divergenceSolve();
    void divergenceSolveIteration(unsigned int fluidModelIndex, Real &avg_density_err);
    void computeDensityAdv(unsigned int fluidModelIndex, unsigned int index, int numParticles, Real h, Real density0);
    void computeDensityChange(unsigned int fluidModelIndex, unsigned int index, Real h);

#ifdef USE_WARMSTART_V
    void warmstartDivergenceSolve(unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
    void warmstartPressureSolve(unsigned int fluidModelIndex);
#endif

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();
    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

    void initParameters() override;

public:
    static int SOLVER_ITERATIONS_V;
    static int MAX_ITERATIONS_V;
    static int MAX_ERROR_V;
    static int USE_DIVERGENCE_SOLVER;

    TimeStepDFSPH();
    ~TimeStepDFSPH() override;

    void step() override;
    void reset() override;

    void resize() override;
};
}  // namespace vox