//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/scene_loader.h"

namespace vox {
/** \brief Class to store the scene configuration that is imported from the scene file.
 */
class SceneConfiguration {
private:
    static SceneConfiguration* m_current;

protected:
    utility::SceneLoader::Scene m_scene;
    std::string m_sceneFile;

public:
    SceneConfiguration();
    SceneConfiguration(const SceneConfiguration&) = delete;
    SceneConfiguration& operator=(const SceneConfiguration&) = delete;
    ~SceneConfiguration();

    // Singleton
    static SceneConfiguration* getCurrent();
    static void setCurrent(SceneConfiguration* sc);
    static bool hasCurrent();

    void setSceneFile(const std::string& file) { m_sceneFile = file; }
    [[nodiscard]] const std::string& getSceneFile() const { return m_sceneFile; }

    utility::SceneLoader::Scene& getScene() { return m_scene; }
};
}  // namespace vox