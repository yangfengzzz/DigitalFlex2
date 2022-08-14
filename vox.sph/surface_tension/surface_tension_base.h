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
/** \brief Base class for all surface tension methods.
 */
class SurfaceTensionBase : public NonPressureForceBase {
protected:
    Real m_surfaceTension;
    Real m_surfaceTensionBoundary;

    void initParameters() override;

public:
    static int SURFACE_TENSION;
    static int SURFACE_TENSION_BOUNDARY;

    explicit SurfaceTensionBase(FluidModel *model);
    ~SurfaceTensionBase() override;
};
}  // namespace vox
