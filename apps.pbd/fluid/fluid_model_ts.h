//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "apps.pbd/fluid/fluid_model.h"

namespace vox {
class TimeStepFluidModel {
protected:
    unsigned int m_velocityUpdateMethod{};

    void clearAccelerations(FluidModel &model);
    void computeXSPHViscosity(FluidModel &model);
    void computeDensities(FluidModel &model);
    void updateTimeStepSizeCFL(FluidModel &model, Real minTimeStepSize, Real maxTimeStepSize);
    void constraintProjection(FluidModel &model);

public:
    TimeStepFluidModel();
    virtual ~TimeStepFluidModel();

    void step(FluidModel &model);
    void reset();

    [[nodiscard]] unsigned int getVelocityUpdateMethod() const { return m_velocityUpdateMethod; }
    void setVelocityUpdateMethod(unsigned int val) { m_velocityUpdateMethod = val; }
};
}  // namespace vox
