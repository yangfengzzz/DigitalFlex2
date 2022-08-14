//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/surface_tension/surface_tension_base.h"

namespace vox {
/** \brief This class implements the surface tension method introduced
 * by Becker and Teschner [BT07].
 *
 * References:
 * - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM
 * SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland,
 * 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
 */
class SurfaceTension_Becker2007 : public SurfaceTensionBase {
public:
    explicit SurfaceTension_Becker2007(FluidModel* model);
    ~SurfaceTension_Becker2007() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Becker2007(model); }

    void step() override;
    void reset() override;
};
}  // namespace vox
