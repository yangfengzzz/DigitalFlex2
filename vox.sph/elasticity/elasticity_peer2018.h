//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/matrix_free_solver.h"
#include "vox.sph/elasticity/elasticity_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief This class implements the implicit SPH formulation for
 * incompressible linearly elastic solids introduced
 * by Peer et al. [PGBT17].
 *
 * References:
 * - [PGBT17] Andreas Peer, Christoph Gissler, Stefan Band, and Matthias Teschner. An implicit SPH formulation for
 * incompressible linearly elastic solids. Computer Graphics Forum, 2017. URL: http://dx.doi.org/10.1111/cgf.13317
 */
class Elasticity_Peer2018 : public ElasticityBase {
protected:
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner>
            Solver;

    // initial particle indices, used to access their original positions
    std::vector<unsigned int> m_current_to_initial_index;
    std::vector<unsigned int> m_initial_to_current_index;
    // initial particle neighborhood
    std::vector<std::vector<unsigned int>> m_initialNeighbors;
    // volumes in rest configuration
    std::vector<Real> m_restVolumes;
    std::vector<Matrix3r> m_rotations;
    std::vector<Vector6r> m_stress;
    std::vector<Matrix3r> m_L;
    std::vector<Matrix3r> m_RL;
    std::vector<Matrix3r> m_F;
    unsigned int m_iterations;
    unsigned int m_maxIter;
    Real m_maxError;
    Real m_alpha;
    Solver m_solver;

    void initValues();
    void computeMatrixL();
    void computeRotations();
    void computeRHS(VectorXr &rhs);

    void initParameters() override;
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    void deferredInit() override;

    //////////////////////////////////////////////////////////////////////////
    // multiplication of symmetric matrix, represented by a 6D vector, and a
    // 3D vector
    //////////////////////////////////////////////////////////////////////////
    static FORCE_INLINE void symMatTimesVec(const Vector6r &M, const Vector3r &v, Vector3r &res) {
        res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
        res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
        res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
    }

public:
    static int ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;
    static int ALPHA;

    explicit Elasticity_Peer2018(FluidModel *model);
    ~Elasticity_Peer2018() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new Elasticity_Peer2018(model); }

    void step() override;
    void reset() override;
    void performNeighborhoodSearchSort() override;

    void saveState(BinaryFileWriter &binWriter) override;
    void loadState(BinaryFileReader &binReader) override;

    static void matrixVecProd(const Real *vec, Real *result, void *userData);
};
}  // namespace vox
