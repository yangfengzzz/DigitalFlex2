//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/non_pressure_force_base.h"

using namespace vox;

NonPressureForceBase::NonPressureForceBase(FluidModel *model) { m_model = model; }

NonPressureForceBase::~NonPressureForceBase() = default;

void NonPressureForceBase::init() { initParameters(); }
