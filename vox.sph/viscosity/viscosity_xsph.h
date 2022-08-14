//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/viscosity/viscosity_base.h"

namespace vox {
/** \brief This class implements the XSPH method descibed by
 * Schechter and Bridson [SB12].
 *
 * References:
 * - [SB12] Hagit Schechter and Robert Bridson. Ghost sph for animating water. ACM Trans. Graph., 31(4):61:1-61:8, July
 * 2012. URL: http://doi.acm.org/10.1145/2185520.2185557
 */
class Viscosity_XSPH : public ViscosityBase {
protected:
    Real m_boundaryViscosity;

    void initParameters() override;

public:
    static int VISCOSITY_COEFFICIENT_BOUNDARY;

    explicit Viscosity_XSPH(FluidModel* model);
    ~Viscosity_XSPH() override;

    void step() override;
    void reset() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_XSPH(model); }
};
}  // namespace vox
