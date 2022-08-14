//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/non_pressure_force_base.h"

namespace vox {
/** \brief Base class for all elasticity methods.
 */
class ElasticityBase : public NonPressureForceBase {
protected:
    Real m_youngsModulus;
    Real m_poissonRatio;
    Vector3r m_fixedBoxMin;
    Vector3r m_fixedBoxMax;

    void initParameters() override;
    void determineFixedParticles();

public:
    static int YOUNGS_MODULUS;
    static int POISSON_RATIO;
    static int FIXED_BOX_MIN;
    static int FIXED_BOX_MAX;

    explicit ElasticityBase(FluidModel *model);
    ~ElasticityBase() override;
};
}  // namespace vox