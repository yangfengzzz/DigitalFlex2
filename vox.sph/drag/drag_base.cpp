//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/drag/drag_base.h"

using namespace vox;

int DragBase::DRAG_COEFFICIENT = -1;

DragBase::DragBase(FluidModel* model) : NonPressureForceBase(model) { m_dragCoefficient = static_cast<Real>(0.01); }

DragBase::~DragBase() = default;

void DragBase::initParameters() {
    NonPressureForceBase::initParameters();

    DRAG_COEFFICIENT = createNumericParameter("drag", "Drag coefficient", &m_dragCoefficient);
    setGroup(DRAG_COEFFICIENT, "Drag force");
    setDescription(DRAG_COEFFICIENT, "Coefficient for the drag force computation");
    auto* rparam = static_cast<RealParameter*>(getParameter(DRAG_COEFFICIENT));
    rparam->setMinValue(0.0);
}
