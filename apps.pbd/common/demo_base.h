//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <AntTweakBar.h>

#include "apps.pbd/common/shader.h"
#include "vox.base/common.h"
#include "vox.base/reflect/parameter_object.h"
#include "vox.pbd/scene_loader.h"
#include "vox.pbd/simulation_model.h"
#include "vox.pbd/time_step.h"

namespace vox {
class DemoBase : public ParameterObject {
protected:
    unsigned int m_numberOfStepsPerRenderUpdate;
    std::string m_exePath;
    std::string m_outputPath;
    std::string m_sceneFile;
    std::string m_sceneName;
    bool m_useCache;
    GLint m_context_major_version{};
    GLint m_context_minor_version{};
    Shader m_shader;
    Shader m_shaderFlat;
    Shader m_shaderTex;
    bool m_doPause;
    Real m_pauseAt;
    bool m_renderTets;
    bool m_renderRefTets;
    bool m_renderContacts;
    bool m_renderAABB;
    bool m_renderSDF;
    int m_renderBVHDepth;
    int m_renderBVHDepthTets;
    Vector3r m_oldMousePos;
    std::vector<unsigned int> m_selectedBodies;
    std::vector<unsigned int> m_selectedParticles;
    utility::SceneLoader *m_sceneLoader;
    utility::SceneLoader::SceneData m_scene;
    float m_jointColor[4] = {0.0f, 0.6f, 0.2f, 1};

    void initParameters() override;

    void initShaders();

    static void selection(const Vector2i &start, const Vector2i &end, void *clientData);
    static void mouseMove(int x, int y, void *clientData);

    void renderTriangleModels();
    void renderTetModels();
    void renderAABB(AABB &aabb);
    void renderSDF(CollisionDetection::CollisionObject *co);
    void renderBallJoint(BallJoint &bj);
    void renderRigidBodyParticleBallJoint(RigidBodyParticleBallJoint &bj);
    void renderBallOnLineJoint(BallOnLineJoint &bj);
    void renderHingeJoint(HingeJoint &hj);
    void renderUniversalJoint(UniversalJoint &uj);
    void renderSliderJoint(SliderJoint &joint);
    void renderTargetPositionMotorSliderJoint(TargetPositionMotorSliderJoint &joint);
    void renderTargetVelocityMotorSliderJoint(TargetVelocityMotorSliderJoint &joint);
    void renderTargetAngleMotorHingeJoint(TargetAngleMotorHingeJoint &hj);
    void renderTargetVelocityMotorHingeJoint(TargetVelocityMotorHingeJoint &hj);
    void renderRigidBodyContact(RigidBodyContactConstraint &cc);
    void renderParticleRigidBodyContact(ParticleRigidBodyContactConstraint &cc);
    void renderSpring(RigidBodySpring &s);
    void renderDistanceJoint(DistanceJoint &j);
    void renderDamperJoint(DamperJoint &j);

public:
    static int PAUSE;
    static int PAUSE_AT;
    static int NUM_STEPS_PER_RENDER;
    static int RENDER_TETS;
    static int RENDER_TETS0;
    static int RENDER_CONTACTS;
    static int RENDER_AABB;
    static int RENDER_SDF;
    static int RENDER_BVH;
    static int RENDER_BVH_TETS;

    DemoBase();
    ~DemoBase() override;

    void init(int argc, char **argv, const char *demoName);

    void createParameterGUI();

    void render();
    void cleanup();

    void readParameters();
    void readScene();
    void reset();

    utility::SceneLoader *getSceneLoader() { return m_sceneLoader; }
    void setSceneLoader(utility::SceneLoader *sceneLoader) { m_sceneLoader = sceneLoader; }

    [[nodiscard]] const std::string &getExePath() const { return m_exePath; }
    [[nodiscard]] const std::string &getSceneFile() const { return m_sceneFile; }
    [[nodiscard]] const std::string &getSceneName() const { return m_sceneName; }

    [[nodiscard]] GLint getContextMajorVersion() const { return m_context_major_version; }
    [[nodiscard]] GLint getContextMinorVersion() const { return m_context_minor_version; }
    Shader &getShader() { return m_shader; }
    Shader &getShaderTex() { return m_shaderTex; }
    Shader &getShaderFlat() { return m_shaderFlat; }
    void shaderTexBegin(const float *col);
    void shaderTexEnd();
    void shaderBegin(const float *col);
    void shaderEnd();
    void shaderFlatBegin(const float *col);
    void shaderFlatEnd();

    std::vector<unsigned int> &getSelectedParticles() { return m_selectedParticles; }
    std::vector<unsigned int> &getSelectedRigidBodies() { return m_selectedBodies; }
    [[nodiscard]] bool getUseCache() const { return m_useCache; }
    void setUseCache(bool val) { m_useCache = val; }
    [[nodiscard]] std::string getOutputPath() const { return m_outputPath; }

    utility::SceneLoader::SceneData &getSceneData() { return m_scene; }
};
}  // namespace vox