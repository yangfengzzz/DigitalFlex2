//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/utilities/matrix_free_solver.h"
#include "vox.sph/viscosity/viscosity_base.h"

namespace vox {
/** \brief This class implements the implicit simulation method for
 * viscous fluids introduced
 * by Peer et al. [PICT15].
 *
 * References:
 * - [PICT15] A. Peer, M. Ihmsen, J. Cornelis, and M. Teschner. An Implicit Viscosity Formulation for SPH Fluids. ACM
 * Trans. Graph., 34(4):1-10, 2015. URL: http://doi.acm.org/10.1145/2766925
 */
class Viscosity_Peer2015 : public ViscosityBase {
protected:
    std::vector<Real> m_density;
    std::vector<Matrix3r> m_targetNablaV;
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner1D> Solver;
    Solver m_solver;
    unsigned int m_iterations;
    unsigned int m_maxIter;
    Real m_maxError;

    void initParameters() override;
    void computeDensities();

public:
    static int ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;

    explicit Viscosity_Peer2015(FluidModel *model);
    ~Viscosity_Peer2015() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new Viscosity_Peer2015(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    static void matrixVecProd(const Real *vec, Real *result, void *userData);
    FORCE_INLINE static void diagonalMatrixElement(unsigned int row, Real &result, void *userData);

    FORCE_INLINE const Matrix3r &getTargetNablaV(const unsigned int i) const { return m_targetNablaV[i]; }

    FORCE_INLINE Matrix3r &getTargetNablaV(const unsigned int i) { return m_targetNablaV[i]; }

    FORCE_INLINE void setTargetNablaV(const unsigned int i, const Matrix3r &val) { m_targetNablaV[i] = val; }
};
}  // namespace vox