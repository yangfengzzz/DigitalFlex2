//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/vorticity/vorticity_base.h"

using namespace vox;

int VorticityBase::VORTICITY_COEFFICIENT = -1;

VorticityBase::VorticityBase(FluidModel* model) : NonPressureForceBase(model) {
    m_vorticityCoeff = static_cast<Real>(0.01);
}

VorticityBase::~VorticityBase() = default;

void VorticityBase::initParameters() {
    NonPressureForceBase::initParameters();

    VORTICITY_COEFFICIENT = createNumericParameter("vorticity", "Vorticity transfer coefficient", &m_vorticityCoeff);
    setGroup(VORTICITY_COEFFICIENT, "Vorticity");
    setDescription(VORTICITY_COEFFICIENT, "Coefficient for the vorticity force computation");
    auto* rparam = static_cast<RealParameter*>(getParameter(VORTICITY_COEFFICIENT));
    rparam->setMinValue(0.0);
}
