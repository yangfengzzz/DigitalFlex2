//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/discrete_grid/discrete_grid.h"
#include "vox.base/reflect/parameter_object.h"
#include "vox.sph/boundary_model.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Base class for the simulation methods.
 */
class TimeStep : public ParameterObject {
public:
    static int SOLVER_ITERATIONS;
    static int MIN_ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;

protected:
    unsigned int m_iterations;
    Real m_maxError;
    unsigned int m_minIterations;
    unsigned int m_maxIterations;

    /** Clear accelerations and add gravitation.
     */
    static void clearAccelerations(unsigned int fluidModelIndex);

    void initParameters() override;

    static void approximateNormal(DiscreteGrid *map, const Eigen::Vector3d &x, Eigen::Vector3d &n, unsigned int dim);
    static void computeVolumeAndBoundaryX(unsigned int fluidModelIndex, unsigned int i, const Vector3r &xi);
    static void computeVolumeAndBoundaryX();
    static void computeDensityAndGradient(unsigned int fluidModelIndex, unsigned int i, const Vector3r &xi);
    static void computeDensityAndGradient();

public:
    TimeStep();
    ~TimeStep() override;

    /** Determine densities of all fluid particles.
     */
    static void computeDensities(unsigned int fluidModelIndex);

    virtual void step() = 0;
    virtual void reset();

    virtual void init();
    virtual void resize() = 0;

    virtual void emittedParticles(FluidModel *model, const unsigned int startIndex){};

    virtual void saveState(BinaryFileWriter &binWriter){};
    virtual void loadState(BinaryFileReader &binReader){};

#ifdef USE_PERFORMANCE_OPTIMIZATION
    void precomputeValues();
#endif
};
}  // namespace vox
