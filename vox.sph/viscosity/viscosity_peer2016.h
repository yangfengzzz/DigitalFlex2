//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/matrix_free_solver.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/viscosity/viscosity_base.h"

namespace vox {
/** \brief This class implements the implicit simulation method for
 * viscous fluids introduced
 * by Peer and Teschner [PGBT17].
 *
 * References:
 * - [PGBT17] Andreas Peer, Christoph Gissler, Stefan Band, and Matthias Teschner. An implicit SPH formulation for
 * incompressible linearly elastic solids. Computer Graphics Forum, 2017. URL: http://dx.doi.org/10.1111/cgf.13317
 */
class Viscosity_Peer2016 : public ViscosityBase {
protected:
    std::vector<Real> m_density;
    std::vector<Matrix3r> m_targetNablaV;
    std::vector<Vector3r> m_omega;
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner1D> Solver;
    Solver m_solverV;
    Solver m_solverOmega;
    unsigned int m_iterationsV;
    unsigned int m_iterationsOmega;
    unsigned int m_maxIterV;
    Real m_maxErrorV;
    unsigned int m_maxIterOmega;
    Real m_maxErrorOmega;

    void initParameters() override;
    void computeDensities();

public:
    static int ITERATIONS_V;
    static int ITERATIONS_OMEGA;
    static int MAX_ITERATIONS_V;
    static int MAX_ERROR_V;
    static int MAX_ITERATIONS_OMEGA;
    static int MAX_ERROR_OMEGA;

    explicit Viscosity_Peer2016(FluidModel *model);
    ~Viscosity_Peer2016() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new Viscosity_Peer2016(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    static void matrixVecProdV(const Real *vec, Real *result, void *userData);
    FORCE_INLINE static void diagonalMatrixElementV(unsigned int row, Real &result, void *userData);

    static void matrixVecProdOmega(const Real *vec, Real *result, void *userData);
    FORCE_INLINE static void diagonalMatrixElementOmega(unsigned int row, Real &result, void *userData);

    FORCE_INLINE const Matrix3r &getTargetNablaV(const unsigned int i) const { return m_targetNablaV[i]; }

    FORCE_INLINE Matrix3r &getTargetNablaV(const unsigned int i) { return m_targetNablaV[i]; }

    FORCE_INLINE void setTargetNablaV(const unsigned int i, const Matrix3r &val) { m_targetNablaV[i] = val; }

    FORCE_INLINE const Vector3r &getOmega(const unsigned int i) const { return m_omega[i]; }

    FORCE_INLINE Vector3r &getOmega(const unsigned int i) { return m_omega[i]; }

    FORCE_INLINE void setOmega(const unsigned int i, const Vector3r &val) { m_omega[i] = val; }
};
}  // namespace vox
