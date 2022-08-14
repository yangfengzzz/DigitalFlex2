//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/scene_configuration.h"

using namespace vox;
using namespace std;

SceneConfiguration* SceneConfiguration::m_current = nullptr;

SceneConfiguration::SceneConfiguration() { m_sceneFile = ""; }

SceneConfiguration::~SceneConfiguration() {
    for (auto& boundaryModel : m_scene.boundaryModels) delete boundaryModel;
    m_scene.boundaryModels.clear();

    for (auto& fluidModel : m_scene.fluidModels) delete fluidModel;
    m_scene.fluidModels.clear();

    for (auto& fluidBlock : m_scene.fluidBlocks) delete fluidBlock;
    m_scene.fluidBlocks.clear();

    for (auto& material : m_scene.materials) delete material;
    m_scene.materials.clear();

    for (auto& emitter : m_scene.emitters) delete emitter;
    m_scene.emitters.clear();

    for (auto& animatedField : m_scene.animatedFields) delete animatedField;
    m_scene.animatedFields.clear();

    m_current = nullptr;
}

SceneConfiguration* SceneConfiguration::getCurrent() {
    if (m_current == nullptr) {
        m_current = new SceneConfiguration();
    }
    return m_current;
}

void SceneConfiguration::setCurrent(SceneConfiguration* sc) { m_current = sc; }

bool SceneConfiguration::hasCurrent() { return (m_current != nullptr); }
