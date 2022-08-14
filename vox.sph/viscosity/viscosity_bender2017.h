//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/viscosity/viscosity_base.h"

namespace vox {
/** \brief This class implements the implicit simulation method for
 * viscous fluids introduced
 * by Bender and Koschier [BK17].
 *
 * References:
 * - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on
 * Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
 */
class Viscosity_Bender2017 : public ViscosityBase {
protected:
    std::vector<Vector6r> m_targetStrainRate;
    std::vector<Matrix6r> m_viscosityFactor;
    std::vector<Vector6r> m_viscosityLambda;
    unsigned int m_iterations;
    unsigned int m_maxIter;
    Real m_maxError;

    void initParameters() override;

public:
    static int ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;

    explicit Viscosity_Bender2017(FluidModel* model);
    ~Viscosity_Bender2017() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_Bender2017(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    void computeTargetStrainRate();
    void computeViscosityFactor();

    /** Matrix product
     */
    static FORCE_INLINE void viscoGradientMultTransposeRightOpt(const Eigen::Matrix<Real, 6, 3>& a,
                                                                const Eigen::Matrix<Real, 6, 3>& b,
                                                                Matrix6r& c) {
        // a(0,1), a(0,2), a(1,0), a(1,2), a(2,0), a(2,1), a(3,2), a(4, 1), a(5, 0) = 0
        c(0, 0) = a(0, 0) * b(0, 0);
        c(0, 1) = 0.0;
        c(0, 2) = 0.0;
        c(0, 3) = a(0, 0) * b(3, 0);
        c(0, 4) = a(0, 0) * b(4, 0);
        c(0, 5) = 0.0;

        c(1, 0) = 0.0;
        c(1, 1) = a(1, 1) * b(1, 1);
        c(1, 2) = 0.0;
        c(1, 3) = a(1, 1) * b(3, 1);
        c(1, 4) = 0.0;
        c(1, 5) = a(1, 1) * b(5, 1);

        c(2, 0) = 0.0;
        c(2, 1) = 0.0;
        c(2, 2) = a(2, 2) * b(2, 2);
        c(2, 3) = 0.0;
        c(2, 4) = a(2, 2) * b(4, 2);
        c(2, 5) = a(2, 2) * b(5, 2);

        c(3, 0) = a(3, 0) * b(0, 0);
        c(3, 1) = a(3, 1) * b(1, 1);
        c(3, 2) = 0.0;
        c(3, 3) = a(3, 0) * b(3, 0) + a(3, 1) * b(3, 1);
        c(3, 4) = a(3, 0) * b(4, 0);
        c(3, 5) = a(3, 1) * b(5, 1);

        c(4, 0) = a(4, 0) * b(0, 0);
        c(4, 1) = 0.0;
        c(4, 2) = a(4, 2) * b(2, 2);
        c(4, 3) = a(4, 0) * b(3, 0);
        c(4, 4) = a(4, 0) * b(4, 0) + a(4, 2) * b(4, 2);
        c(4, 5) = a(4, 2) * b(5, 2);

        c(5, 0) = 0.0;
        c(5, 1) = a(5, 1) * b(1, 1);
        c(5, 2) = a(5, 2) * b(2, 2);
        c(5, 3) = a(5, 1) * b(3, 1);
        c(5, 4) = a(5, 2) * b(4, 2);
        c(5, 5) = a(5, 1) * b(5, 1) + a(5, 2) * b(5, 2);
    }

    [[nodiscard]] FORCE_INLINE const Vector6r& getTargetStrainRate(const unsigned int i) const {
        return m_targetStrainRate[i];
    }

    FORCE_INLINE Vector6r& getTargetStrainRate(const unsigned int i) { return m_targetStrainRate[i]; }

    FORCE_INLINE void setTargetStrainRate(const unsigned int i, const Vector6r& val) { m_targetStrainRate[i] = val; }

    [[nodiscard]] FORCE_INLINE const Matrix6r& getViscosityFactor(const unsigned int i) const {
        return m_viscosityFactor[i];
    }

    FORCE_INLINE Matrix6r& getViscosityFactor(const unsigned int i) { return m_viscosityFactor[i]; }

    FORCE_INLINE void setViscosityFactor(const unsigned int i, const Matrix6r& val) { m_viscosityFactor[i] = val; }

    [[nodiscard]] FORCE_INLINE const Vector6r& getViscosityLambda(const unsigned int i) const {
        return m_viscosityLambda[i];
    }

    FORCE_INLINE Vector6r& getViscosityLambda(const unsigned int i) { return m_viscosityLambda[i]; }

    FORCE_INLINE void setViscosityLambda(const unsigned int i, const Vector6r& val) { m_viscosityLambda[i] = val; }
};
}  // namespace vox
