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
 * by Kugelstadt et al. [KBF+21].
 *
 * References:
 * - [KBF+21] Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys Jeske, Fabian Löschner,
 * Andreas Longva. Fast Corotated Elastic SPH Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM on
 * Computer Graphics and Interactive Techniques, 2021. URL: http://dx.doi.org/10.1145/3480142
 */
class Elasticity_Kugelstadt2021 : public ElasticityBase {
protected:
    struct Factorization {
        Real m_dt;
        Real m_mu;
        Eigen::SparseMatrix<Real, Eigen::RowMajor> m_DT_K;
        Eigen::SparseMatrix<Real, Eigen::RowMajor> m_D;
        Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matHTH;
        Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matL;
        Eigen::SparseMatrix<Real, Eigen::ColMajor> m_matLT;
        Eigen::VectorXi m_permInd;
        Eigen::VectorXi m_permInvInd;
    };

    struct ElasticObject {
        std::string m_md5;
        std::vector<unsigned int> m_particleIndices;
        unsigned int m_nFixed{};

        std::shared_ptr<Factorization> m_factorization;
        std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_f;
        std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_sol;
        std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_RHS;
        std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_RHS_perm;
        std::vector<Quaternionr, Eigen::aligned_allocator<Quaternionr>> m_quats;

        ElasticObject() { m_factorization = nullptr; }
        ~ElasticObject() { m_factorization = nullptr; }
    };

    // initial particle indices, used to access their original positions
    std::vector<unsigned int> m_current_to_initial_index;
    std::vector<unsigned int> m_initial_to_current_index;
    // initial particle neighborhood
    std::vector<std::vector<unsigned int>> m_initialNeighbors;
    // volumes in rest configuration
    std::vector<Real> m_restVolumes;
    std::vector<Matrix3r> m_rotations;
    std::vector<Real> m_stress;
    std::vector<Matrix3r> m_L;
    std::vector<Matrix3r> m_F;
    std::vector<Vector3r> m_vDiff;
    std::vector<Matrix3r> m_RL;
    unsigned int m_iterationsV;
    unsigned int m_maxIterV;
    Real m_maxErrorV;
    Real m_alpha;
    int m_maxNeighbors;
    unsigned int m_totalNeighbors;
    std::vector<ElasticObject*> m_objects;
    Real m_lambda;
    Real m_mu;

    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_precomp_RL_gradW;
    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_precomp_L_gradW;
    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> m_precomp_RLj_gradW;
    std::vector<unsigned int> m_precomputed_indices;

    typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> SolverLLT;
    typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> Solver;

    Solver m_solver;
    void computeRHS(VectorXr& rhs);

    std::string computeMD5(unsigned int objIndex);
    void initValues();
    void initSystem();
    void initFactorization(std::shared_ptr<Factorization> factorization,
                           std::vector<unsigned int>& particleIndices,
                           unsigned int nFixed,
                           Real dt,
                           Real mu);
    void findObjects();
    void computeMatrixL();
    void precomputeValues();

    void stepElasticitySolver();
    void stepVolumeSolver();

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
    static FORCE_INLINE void symMatTimesVec(const Vector6r& M, const Vector3r& v, Vector3r& res) {
        res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
        res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
        res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
    }

    void rotationMatricesToAVXQuaternions();

public:
    static int ITERATIONS_V;
    static int MAX_ITERATIONS_V;
    static int MAX_ERROR_V;
    static int ALPHA;
    static int MAX_NEIGHBORS;

    explicit Elasticity_Kugelstadt2021(FluidModel* model);
    ~Elasticity_Kugelstadt2021() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new Elasticity_Kugelstadt2021(model); }

    void step() override;
    void reset() override;
    void performNeighborhoodSearchSort() override;

    void saveState(BinaryFileWriter& binWriter) override;
    void loadState(BinaryFileReader& binReader) override;

    static void matrixVecProd(const Real* vec, Real* result, void* userData);

    void computeRotations();
};
}  // namespace vox
