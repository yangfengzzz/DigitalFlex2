//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/drag/drag_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief This class implements the drag force computation introduced
 * by Macklin et al. [MMCK14].
 *
 * References:
 * - [MMCK14] Miles Macklin, Matthias Müller, Nuttapong Chentanez, and Tae-Yong Kim. Unified Particle Physics for
 * Real-Time Applications. ACM Trans. Graph., 33(4):1-12, 2014. URL: http://doi.acm.org/10.1145/2601097.2601152
 */
class DragForce_Macklin2014 : public DragBase {
public:
    explicit DragForce_Macklin2014(FluidModel* model);
    ~DragForce_Macklin2014() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new DragForce_Macklin2014(model); }

    void step() override;
    void reset() override;
};
}  // namespace vox