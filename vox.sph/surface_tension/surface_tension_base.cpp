//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/surface_tension/surface_tension_base.h"

using namespace vox;

int SurfaceTensionBase::SURFACE_TENSION = -1;
int SurfaceTensionBase::SURFACE_TENSION_BOUNDARY = -1;

SurfaceTensionBase::SurfaceTensionBase(FluidModel* model) : NonPressureForceBase(model) {
    m_surfaceTension = 0.05;
    m_surfaceTensionBoundary = 0.01;
}

SurfaceTensionBase::~SurfaceTensionBase() = default;

void SurfaceTensionBase::initParameters() {
    NonPressureForceBase::initParameters();

    SURFACE_TENSION = createNumericParameter("surfaceTension", "Surface tension coefficient", &m_surfaceTension);
    setGroup(SURFACE_TENSION, "Surface tension");
    setDescription(SURFACE_TENSION, "Coefficient for the surface tension computation");
    auto* rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION));
    rparam->setMinValue(0.0);

    SURFACE_TENSION_BOUNDARY = createNumericParameter("surfaceTensionBoundary", "Boundary surface tension coefficient",
                                                      &m_surfaceTensionBoundary);
    setGroup(SURFACE_TENSION_BOUNDARY, "Surface tension");
    setDescription(SURFACE_TENSION_BOUNDARY, "Coefficient for the surface tension computation at the boundary");
    rparam = static_cast<RealParameter*>(getParameter(SURFACE_TENSION_BOUNDARY));
    rparam->setMinValue(0.0);
}
