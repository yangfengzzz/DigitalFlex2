//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "apps.pbd/elastic_rods/position_based_elastic_rods_model.h"
#include "vox.pbd/time_step_controller.h"

namespace vox {
class PositionBasedElasticRodsTSC : public TimeStepController {
protected:
    Real m_damping;
    virtual void clearAccelerations(SimulationModel &model);

public:
    PositionBasedElasticRodsTSC();
    ~PositionBasedElasticRodsTSC() override;

    void step(SimulationModel &model) override;

    [[nodiscard]] Real getDamping() const { return m_damping; }
    void setDamping(Real val) { m_damping = val; }
};
}  // namespace vox