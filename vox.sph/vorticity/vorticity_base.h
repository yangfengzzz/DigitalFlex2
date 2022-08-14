//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/non_pressure_force_base.h"

namespace vox {
/** \brief Base class for all vorticity methods.
 */
class VorticityBase : public NonPressureForceBase {
protected:
    Real m_vorticityCoeff;

    void initParameters() override;

public:
    static int VORTICITY_COEFFICIENT;

    explicit VorticityBase(FluidModel *model);
    ~VorticityBase() override;
};
}  // namespace vox
