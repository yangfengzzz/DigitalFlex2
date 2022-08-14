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
 * by Gissler et al. [GBP+17].
 *
 * References:
 * - [GPB+17] Christoph Gissler, Stefan Band, Andreas Peer, Markus Ihmsen, and Matthias Teschner. Approximate air-fluid
 * interactions for SPH. In Virtual Reality Interactions and Physical Simulations, 1-10. April 2017. URL:
 * http://dx.doi.org/10.2312/vriphys.20171081
 */
class DragForce_Gissler2017 : public DragBase {
protected:
    const Real rho_a = static_cast<Real>(1.2041);
    const Real sigma = static_cast<Real>(0.0724);
    const Real mu_l = static_cast<Real>(0.00102);
    const Real C_F = static_cast<Real>(1.0 / 3.0);
    const Real C_k = static_cast<Real>(8.0);
    const Real C_d = static_cast<Real>(5.0);
    const Real C_b = static_cast<Real>(0.5);
    const Real mu_a = static_cast<Real>(0.00001845);

public:
    explicit DragForce_Gissler2017(FluidModel* model);
    ~DragForce_Gissler2017() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new DragForce_Gissler2017(model); }

    void step() override;
    void reset() override;
};
}  // namespace vox