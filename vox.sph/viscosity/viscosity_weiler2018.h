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

#define USE_BLOCKDIAGONAL_PRECONDITIONER

namespace vox {
/** \brief This class implements the implicit Laplace viscosity method introduced
 * by Weiler et al. 2018 [WKBB18].
 *
 * References:
 * - [WKBB18] Marcel Weiler, Dan Koschier, Magnus Brand, and Jan Bender. A physically consistent implicit viscosity
 * solver for SPH fluids. Computer Graphics Forum (Eurographics), 2018. URL: https://doi.org/10.1111/cgf.13349
 */
class Viscosity_Weiler2018 : public ViscosityBase {
protected:
    Real m_boundaryViscosity;
    unsigned int m_maxIter;
    Real m_maxError;
    unsigned int m_iterations;
    std::vector<Vector3r> m_vDiff;
    Real m_tangentialDistanceFactor;

#ifdef USE_BLOCKDIAGONAL_PRECONDITIONER
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, BlockJacobiPreconditioner3D>
            Solver;
    FORCE_INLINE static void diagonalMatrixElement(unsigned int row, Matrix3r &result, void *userData);
#else
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner3D> Solver;
    FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData);
#endif
    Solver m_solver;

    void initParameters() override;

public:
    static int ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;
    static int VISCOSITY_COEFFICIENT_BOUNDARY;

    explicit Viscosity_Weiler2018(FluidModel *model);
    ~Viscosity_Weiler2018() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new Viscosity_Weiler2018(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    static void matrixVecProd(const Real *vec, Real *result, void *userData);

    FORCE_INLINE const Vector3r &getVDiff(const unsigned int i) const { return m_vDiff[i]; }

    FORCE_INLINE Vector3r &getVDiff(const unsigned int i) { return m_vDiff[i]; }

    FORCE_INLINE void setVDiff(const unsigned int i, const Vector3r &val) { m_vDiff[i] = val; }

    void computeRHS(VectorXr &b, VectorXr &g);

    void applyForces(const VectorXr &x);
};
}  // namespace vox
