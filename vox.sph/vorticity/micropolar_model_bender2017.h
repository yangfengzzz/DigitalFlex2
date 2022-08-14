//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/vorticity/vorticity_base.h"

namespace vox {
/** \brief This class implements the micropolar material model introduced
 * by Bender et al. [BKKW17].
 *
 * References:
 * - [BKKW17] Jan Bender, Dan Koschier, Tassilo Kugelstadt, and Marcel Weiler. A micropolar material model for turbulent
 * SPH fluids. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '17. ACM, 2017. URL:
 * http://doi.acm.org/10.1145/3099564.3099578
 */
class MicropolarModel_Bender2017 : public VorticityBase {
protected:
    std::vector<Vector3r> m_angularAcceleration;
    std::vector<Vector3r> m_omega;
    Real m_viscosityOmega;
    Real m_inertiaInverse;

    void initParameters() override;

public:
    static int VISCOSITY_OMEGA;
    static int INERTIA_INVERSE;

    explicit MicropolarModel_Bender2017(FluidModel* model);
    ~MicropolarModel_Bender2017() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new MicropolarModel_Bender2017(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    [[nodiscard]] FORCE_INLINE const Vector3r& getAngularAcceleration(const unsigned int i) const {
        return m_angularAcceleration[i];
    }

    FORCE_INLINE Vector3r& getAngularAcceleration(const unsigned int i) { return m_angularAcceleration[i]; }

    FORCE_INLINE void setAngularAcceleration(const unsigned int i, const Vector3r& val) {
        m_angularAcceleration[i] = val;
    }

    [[nodiscard]] FORCE_INLINE const Vector3r& getAngularVelocity(const unsigned int i) const { return m_omega[i]; }

    FORCE_INLINE Vector3r& getAngularVelocity(const unsigned int i) { return m_omega[i]; }

    FORCE_INLINE void setAngularVelocity(const unsigned int i, const Vector3r& val) { m_omega[i] = val; }
};
}  // namespace vox