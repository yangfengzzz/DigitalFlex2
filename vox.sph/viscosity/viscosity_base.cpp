//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/viscosity/viscosity_base.h"

using namespace vox;

int ViscosityBase::VISCOSITY_COEFFICIENT = -1;

ViscosityBase::ViscosityBase(FluidModel* model) : NonPressureForceBase(model) { m_viscosity = 0.01; }

ViscosityBase::~ViscosityBase() = default;

void ViscosityBase::initParameters() {
    NonPressureForceBase::initParameters();

    VISCOSITY_COEFFICIENT = createNumericParameter("viscosity", "Viscosity coefficient", &m_viscosity);
    setGroup(VISCOSITY_COEFFICIENT, "Viscosity");
    setDescription(VISCOSITY_COEFFICIENT, "Coefficient for the viscosity force computation");
    auto* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT));
    rparam->setMinValue(0.0);
}
