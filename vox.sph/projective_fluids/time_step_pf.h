//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/matrix_free_solver.h"
#include "vox.base/sph_kernels.h"
#include "vox.sph/projective_fluids/simulation_data_pf.h"
#include "vox.sph/time_step.h"

// since all diagonal blocks are 3x3 diagonal matrices, a diagonal preconditioner does suffice
#define PD_USE_DIAGONAL_PRECONDITIONER

namespace vox {
/** \brief This class implements the Projective Fluids approach introduced
 * by Weiler, Koschier and Bender [WKB16].
 *
 * References:
 * - [WKB16] Marcel Weiler, Dan Koschier, and Jan Bender. Projective fluids. In Proceedings of the 9th International
 * Conference on Motion in Games, MIG '16, 79-84. New York, NY, USA, 2016. ACM. URL:
 * http://doi.acm.org/10.1145/2994258.2994282
 */
class TimeStepPF : public TimeStep {
protected:
    using VectorXr = Eigen::Matrix<Real, -1, 1>;
    using VectorXrMap = Eigen::Map<VectorXr>;

#ifdef PD_USE_DIAGONAL_PRECONDITIONER
    using Solver = Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner3D>;
    FORCE_INLINE static void diagonalMatrixElement(unsigned int row, Vector3r &result, void *userData);
    void preparePreconditioner();
#else
    using Solver =
            Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner>;
#endif

    SimulationDataPF m_simulationData;
    Solver m_solver;
    Real m_stiffness;
    unsigned int m_counter;
    unsigned int m_numActiveParticlesTotal;

    void initialGuessForPositions(unsigned int fluidModelIndex);
    void solvePDConstraints();
    void updatePositionsAndVelocity(const VectorXr &x);
    static void addAccellerationToVelocity();

    void matrixFreeRHS(const VectorXr &x, VectorXr &result);

    /** Perform the neighborhood search for all fluid particles.
     */
    void performNeighborhoodSearch();
    void emittedParticles(FluidModel *model, unsigned int startIndex) override;

    void initParameters() override;

public:
    static int STIFFNESS;

    TimeStepPF();
    ~TimeStepPF() override;

    void step() override;
    void reset() override;
    void resize() override;

    static void matrixVecProd(const Real *vec, Real *result, void *userData);
};
}  // namespace vox
