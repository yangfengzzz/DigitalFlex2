//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.base/reflect/parameter_object.h"
#include "vox.sph/fluid_model.h"

namespace vox {
class DebugTools : public ParameterObject {
protected:
    bool m_determineThreadIds;
    std::vector<std::vector<unsigned int>> m_threadIds;
    bool m_determineNumNeighbors{};
    std::vector<std::vector<unsigned int>> m_numNeighbors;
    bool m_determineVelocityChanges{};
    std::vector<std::vector<Vector3r>> m_vOld;
    std::vector<std::vector<Vector3r>> m_velocityChanges;

    void initParameters() override;

    void determineThreadIds();
    void determineNumNeighbors();
    void determineVelocityChanges();

public:
    static int DETERMINE_THREAD_IDS;
    static int DETERMINE_NUM_NEIGHBORS;
    static int DETERMINE_VELOCITY_CHANGES;

    DebugTools();
    ~DebugTools() override;

    void init();
    void cleanup();

    void step();
    void reset();

    void performNeighborhoodSearchSort();
    void emittedParticles(FluidModel* model, unsigned int startIndex);
};
}  // namespace vox