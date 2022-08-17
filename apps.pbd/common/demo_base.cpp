//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "apps.pbd/common/demo_base.h"

#include "apps.pbd/common/mini_gl.h"
#include "apps.pbd/common/selection.h"
#include "apps.pbd/common/tweakbar_parameters.h"
#include "apps.pbd/common/visualization.h"
#include "vox.base/file_system.h"
#include "vox.base/logging.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.pbd/distance_field_collision_detection.h"
#include "vox.pbd/simulation.h"

INIT_TIMING

using namespace vox;
using namespace std;

int DemoBase::PAUSE = -1;
int DemoBase::PAUSE_AT = -1;
int DemoBase::NUM_STEPS_PER_RENDER = -1;
int DemoBase::RENDER_TETS = -1;
int DemoBase::RENDER_TETS0 = -1;
int DemoBase::RENDER_CONTACTS = -1;
int DemoBase::RENDER_AABB = -1;
int DemoBase::RENDER_SDF = -1;
int DemoBase::RENDER_BVH = -1;
int DemoBase::RENDER_BVH_TETS = -1;

DemoBase::DemoBase() {
    m_sceneLoader = nullptr;
    m_numberOfStepsPerRenderUpdate = 8;
    m_renderContacts = false;
    m_renderAABB = false;
    m_renderSDF = false;
    m_renderBVHDepth = -1;
    m_renderBVHDepthTets = -1;
    m_renderRefTets = false;
    m_renderTets = false;
    m_sceneFile = "";
    m_sceneName = "";
    m_doPause = true;
    m_pauseAt = -1.0;
    m_useCache = true;
    m_oldMousePos.setZero();
}

DemoBase::~DemoBase() { delete m_sceneLoader; }

void DemoBase::initParameters() {
    ParameterObject::initParameters();

    PAUSE = createBoolParameter("pause", "Pause", &m_doPause);
    setGroup(PAUSE, "Simulation");
    setDescription(PAUSE, "Pause simulation.");
    setHotKey(PAUSE, "space");

    PAUSE_AT = createNumericParameter("pauseAt", "Pause simulation at", &m_pauseAt);
    setGroup(PAUSE_AT, "Simulation");
    setDescription(PAUSE_AT,
                   "Pause simulation at the given time. When the value is negative, the simulation is not paused.");

    NUM_STEPS_PER_RENDER = createNumericParameter("numberOfStepsPerRenderUpdate", "# time steps / update",
                                                  &m_numberOfStepsPerRenderUpdate);
    setGroup(NUM_STEPS_PER_RENDER, "Visualization");
    setDescription(NUM_STEPS_PER_RENDER,
                   "Pause simulation at the given time. When the value is negative, the simulation is not paused.");
    static_cast<NumericParameter<unsigned int> *>(getParameter(NUM_STEPS_PER_RENDER))->setMinValue(1);

    RENDER_TETS = createBoolParameter("renderTets", "Render tet model", &m_renderTets);
    setGroup(RENDER_TETS, "Visualization");
    setDescription(RENDER_TETS, "Render tet model.");

    RENDER_TETS0 = createBoolParameter("renderTets0", "Render ref. tet model", &m_renderRefTets);
    setGroup(RENDER_TETS0, "Visualization");
    setDescription(RENDER_TETS0, "Render ref. tet model.");

    RENDER_CONTACTS = createBoolParameter("renderContacts", "Render contacts", &m_renderContacts);
    setGroup(RENDER_CONTACTS, "Visualization");
    setDescription(RENDER_CONTACTS, "Render contact points and normals.");

    RENDER_AABB = createBoolParameter("renderAABB", "Render AABBs", &m_renderAABB);
    setGroup(RENDER_AABB, "Visualization");
    setDescription(RENDER_AABB, "Render AABBs.");

    RENDER_SDF = createBoolParameter("renderSDF", "Render SDFs", &m_renderSDF);
    setGroup(RENDER_SDF, "Visualization");
    setDescription(RENDER_SDF, "Render SDFs.");

    RENDER_BVH = createNumericParameter("renderBVHDepth", "Render BVH depth", &m_renderBVHDepth);
    setGroup(RENDER_BVH, "Visualization");
    setDescription(RENDER_BVH, "Render BVH until given depth.");
    static_cast<NumericParameter<int> *>(getParameter(RENDER_BVH))->setMinValue(-1);

    RENDER_BVH_TETS = createNumericParameter("renderBVHDepthTets", "Render BVH depth (tets)", &m_renderBVHDepthTets);
    setGroup(RENDER_BVH_TETS, "Visualization");
    setDescription(RENDER_BVH_TETS, "Render BVH (tets) until given depth.");
    static_cast<NumericParameter<int> *>(getParameter(RENDER_BVH_TETS))->setMinValue(-1);
}

void DemoBase::createParameterGUI() {
    TwRemoveAllVars(MiniGL::getTweakBar());
    TweakBarParameters::cleanup();

    MiniGL::initTweakBarParameters();

    TweakBarParameters::createParameterGUI();
    TweakBarParameters::createParameterObjectGUI(this);
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getModel());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getTimeStep());
}

void DemoBase::cleanup() {
    m_scene.m_rigidBodyData.clear();
    m_scene.m_rigidBodyData.clear();
    m_scene.m_triangleModelData.clear();
    m_scene.m_tetModelData.clear();
    m_scene.m_ballJointData.clear();
    m_scene.m_ballOnLineJointData.clear();
    m_scene.m_hingeJointData.clear();
    m_scene.m_universalJointData.clear();
    m_scene.m_sliderJointData.clear();
    m_scene.m_rigidBodyParticleBallJointData.clear();
    m_scene.m_targetAngleMotorHingeJointData.clear();
    m_scene.m_targetVelocityMotorHingeJointData.clear();
    m_scene.m_targetPositionMotorSliderJointData.clear();
    m_scene.m_targetVelocityMotorSliderJointData.clear();
    m_scene.m_rigidBodySpringData.clear();
    m_scene.m_distanceJointData.clear();
    m_scene.m_damperJointData.clear();
}

void DemoBase::init(int argc, char **argv, const char *demoName) {
    initParameters();
    m_exePath = utility::FileSystem::getProgramPath();

    m_sceneFile = "";
    setUseCache(true);
    for (int i = 1; i < argc; i++) {
        string argStr = argv[i];
        if (argStr == "--no-cache")
            setUseCache(false);
        else {
            m_sceneFile = string(argv[i]);
            if (utility::FileSystem::isRelativePath(m_sceneFile))
                m_sceneFile = utility::FileSystem::normalizePath(m_exePath + "/" + m_sceneFile);
        }
    }

    m_outputPath = utility::FileSystem::normalizePath(getExePath() + "/output/" +
                                                      utility::FileSystem::getFileName(m_sceneFile));

#ifdef DL_OUTPUT
    std::string sceneFilePath = FileSystem::normalizePath(m_outputPath + "/scene");
    FileSystem::makeDirs(sceneFilePath);
    FileSystem::copyFile(m_sceneFile, sceneFilePath + "/" + FileSystem::getFileNameWithExt(m_sceneFile));

    std::string progFilePath = FileSystem::normalizePath(m_outputPath + "/program");
    FileSystem::makeDirs(progFilePath);
    FileSystem::copyFile(argv[0], progFilePath + "/" + FileSystem::getFileNameWithExt(argv[0]));
#endif

    std::string logPath = utility::FileSystem::normalizePath(m_outputPath + "/log");
    utility::FileSystem::makeDirs(logPath);

    //    LOG_DEBUG << "Git refspec: " << GIT_REFSPEC;
    //    LOG_DEBUG << "Git SHA1:    " << GIT_SHA1;
    //    LOG_DEBUG << "Git status:  " << GIT_LOCAL_STATUS;
    //    LOG_DEBUG << "Host name:   " << SystemInfo::getHostName();
    //    LOG_INFO << "PositionBasedDynamics " << PBD_VERSION;

    // OpenGL
    MiniGL::init(argc, argv, 1280, 1024, demoName);
    MiniGL::initLights();
    MiniGL::initTexture();
    MiniGL::getOpenGLVersion(m_context_major_version, m_context_minor_version);
    MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(0.0, 3.0, 8.0), Vector3r(0.0, 0.0, 0.0));
    MiniGL::setSelectionFunc(selection, this);

    if (MiniGL::checkOpenGLVersion(3, 3)) initShaders();

    MiniGL::addReshapeFunc([](int width, int height) { TwWindowSize(width, height); });
    MiniGL::addKeyboardFunc(
            [](int key, int scancode, int action, int mods) -> bool { return TwEventKeyGLFW(key, action); });
    MiniGL::addCharFunc([](int key, int action) -> bool { return TwEventCharGLFW(key, action); });
    MiniGL::addMousePressFunc(
            [](int button, int action, int mods) -> bool { return TwEventMouseButtonGLFW(button, action); });
    MiniGL::addMouseMoveFunc([](int x, int y) -> bool { return TwEventMousePosGLFW(x, y); });
    MiniGL::addMouseWheelFunc(
            [](int pos, double xoffset, double yoffset) -> bool { return TwEventMouseWheelGLFW(pos); });
}

void DemoBase::readScene() {
    if (m_sceneLoader == nullptr) m_sceneLoader = new utility::SceneLoader();
    if (!m_sceneFile.empty())
        m_sceneLoader->readScene(m_sceneFile, m_scene);
    else
        return;

    m_sceneName = m_scene.m_sceneName;
}

void DemoBase::initShaders() {
    std::string vertFile = m_exePath + "/resources/shaders/vs_smooth.glsl";
    std::string fragFile = m_exePath + "/resources/shaders/fs_smooth.glsl";
    m_shader.compileShaderFile(GL_VERTEX_SHADER, vertFile);
    m_shader.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
    m_shader.createAndLinkProgram();
    m_shader.begin();
    m_shader.addUniform("modelview_matrix");
    m_shader.addUniform("projection_matrix");
    m_shader.addUniform("surface_color");
    m_shader.addUniform("shininess");
    m_shader.addUniform("specular_factor");
    m_shader.end();

    vertFile = m_exePath + "/resources/shaders/vs_smoothTex.glsl";
    fragFile = m_exePath + "/resources/shaders/fs_smoothTex.glsl";
    m_shaderTex.compileShaderFile(GL_VERTEX_SHADER, vertFile);
    m_shaderTex.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
    m_shaderTex.createAndLinkProgram();
    m_shaderTex.begin();
    m_shaderTex.addUniform("modelview_matrix");
    m_shaderTex.addUniform("projection_matrix");
    m_shaderTex.addUniform("surface_color");
    m_shaderTex.addUniform("shininess");
    m_shaderTex.addUniform("specular_factor");
    m_shaderTex.end();

    vertFile = m_exePath + "/resources/shaders/vs_flat.glsl";
    std::string geomFile = m_exePath + "/resources/shaders/gs_flat.glsl";
    fragFile = m_exePath + "/resources/shaders/fs_flat.glsl";
    m_shaderFlat.compileShaderFile(GL_VERTEX_SHADER, vertFile);
    //    m_shaderFlat.compileShaderFile(GL_GEOMETRY_SHADER, geomFile);
    m_shaderFlat.compileShaderFile(GL_FRAGMENT_SHADER, fragFile);
    m_shaderFlat.createAndLinkProgram();
    m_shaderFlat.begin();
    m_shaderFlat.addUniform("modelview_matrix");
    m_shaderFlat.addUniform("projection_matrix");
    m_shaderFlat.addUniform("surface_color");
    m_shaderFlat.addUniform("shininess");
    m_shaderFlat.addUniform("specular_factor");
    m_shaderFlat.end();
}

void DemoBase::shaderTexBegin(const float *col) {
    m_shaderTex.begin();
    glUniform1f(m_shaderTex.getUniform("shininess"), 5.0f);
    glUniform1f(m_shaderTex.getUniform("specular_factor"), 0.2f);

    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    glUniformMatrix4fv(m_shaderTex.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
    GLfloat pmatrix[16];
    glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
    glUniformMatrix4fv(m_shaderTex.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
    glUniform3fv(m_shaderTex.getUniform("surface_color"), 1, col);
}

void DemoBase::shaderTexEnd() { m_shaderTex.end(); }

void DemoBase::shaderBegin(const float *col) {
    m_shader.begin();
    glUniform1f(m_shader.getUniform("shininess"), 5.0f);
    glUniform1f(m_shader.getUniform("specular_factor"), 0.2f);

    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    glUniformMatrix4fv(m_shader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
    GLfloat pmatrix[16];
    glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
    glUniformMatrix4fv(m_shader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
    glUniform3fv(m_shader.getUniform("surface_color"), 1, col);
}

void DemoBase::shaderEnd() { m_shader.end(); }

void DemoBase::shaderFlatBegin(const float *col) {
    m_shaderFlat.begin();
    glUniform1f(m_shaderFlat.getUniform("shininess"), 5.0f);
    glUniform1f(m_shaderFlat.getUniform("specular_factor"), 0.2f);

    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    glUniformMatrix4fv(m_shaderFlat.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
    GLfloat pmatrix[16];
    glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
    glUniformMatrix4fv(m_shaderFlat.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
    glUniform3fv(m_shaderFlat.getUniform("surface_color"), 1, col);
}

void DemoBase::shaderFlatEnd() { m_shaderFlat.end(); }

void DemoBase::readParameters() {
    m_sceneLoader->readParameterObject(this);
    m_sceneLoader->readParameterObject(Simulation::getCurrent());
    m_sceneLoader->readParameterObject(Simulation::getCurrent()->getModel());
    m_sceneLoader->readParameterObject(Simulation::getCurrent()->getTimeStep());
}

void DemoBase::render() {
    float gridColor[4] = {0.2f, 0.2f, 0.2f, 1.0f};
    MiniGL::drawGrid_xz(gridColor);

    MiniGL::coordinateSystem();

    // Draw sim model
    SimulationModel *model = Simulation::getCurrent()->getModel();
    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    SimulationModel::ConstraintVector &constraints = model->getConstraints();
    SimulationModel::RigidBodyContactConstraintVector &rigidBodyContacts = model->getRigidBodyContactConstraints();
    SimulationModel::ParticleRigidBodyContactConstraintVector &particleRigidBodyContacts =
            model->getParticleRigidBodyContactConstraints();

    float selectionColor[4] = {0.8f, 0.0f, 0.0f, 1};
    float surfaceColor[4] = {0.1f, 0.4f, 0.7f, 1};
    float staticColor[4] = {0.5f, 0.5f, 0.5f, 1};

    if (m_renderContacts) {
        for (auto &rigidBodyContact : rigidBodyContacts) renderRigidBodyContact(rigidBodyContact);
        for (auto &particleRigidBodyContact : particleRigidBodyContacts)
            renderParticleRigidBodyContact(particleRigidBodyContact);
    }

    for (size_t i = 0; i < rb.size(); i++) {
        bool selected = false;
        for (unsigned int m_selectedBodie : m_selectedBodies) {
            if (m_selectedBodie == i) selected = true;
        }

        const VertexData &vd = rb[i]->getGeometry().getVertexData();
        const utility::IndexedFaceMesh &mesh = rb[i]->getGeometry().getMesh();

        if (mesh.getFlatShading())
            shaderFlatBegin(staticColor);
        else
            shaderBegin(staticColor);

        if (!selected) {
            if (rb[i]->getMass() == 0.0) {
                glUniform3fv(m_shader.getUniform("surface_color"), 1, staticColor);
                Visualization::drawMesh(vd, mesh, 0, staticColor);
            } else {
                glUniform3fv(m_shader.getUniform("surface_color"), 1, surfaceColor);
                Visualization::drawMesh(vd, mesh, 0, surfaceColor);
            }
        } else {
            glUniform3fv(m_shader.getUniform("surface_color"), 1, selectionColor);
            Visualization::drawMesh(vd, mesh, 0, selectionColor);
        }

        if (mesh.getFlatShading())
            shaderFlatEnd();
        else
            shaderEnd();
    }

    renderTriangleModels();
    renderTetModels();

    for (auto &constraint : constraints) {
        if (constraint->getTypeId() == BallJoint::TYPE_ID) {
            renderBallJoint(*(BallJoint *)constraint);
        } else if (constraint->getTypeId() == BallOnLineJoint::TYPE_ID) {
            renderBallOnLineJoint(*(BallOnLineJoint *)constraint);
        } else if (constraint->getTypeId() == HingeJoint::TYPE_ID) {
            renderHingeJoint(*(HingeJoint *)constraint);
        } else if (constraint->getTypeId() == UniversalJoint::TYPE_ID) {
            renderUniversalJoint(*(UniversalJoint *)constraint);
        } else if (constraint->getTypeId() == SliderJoint::TYPE_ID) {
            renderSliderJoint(*(SliderJoint *)constraint);
        } else if (constraint->getTypeId() == TargetAngleMotorHingeJoint::TYPE_ID) {
            renderTargetAngleMotorHingeJoint(*(TargetAngleMotorHingeJoint *)constraint);
        } else if (constraint->getTypeId() == TargetVelocityMotorHingeJoint::TYPE_ID) {
            renderTargetVelocityMotorHingeJoint(*(TargetVelocityMotorHingeJoint *)constraint);
        } else if (constraint->getTypeId() == TargetPositionMotorSliderJoint::TYPE_ID) {
            renderTargetPositionMotorSliderJoint(*(TargetPositionMotorSliderJoint *)constraint);
        } else if (constraint->getTypeId() == TargetVelocityMotorSliderJoint::TYPE_ID) {
            renderTargetVelocityMotorSliderJoint(*(TargetVelocityMotorSliderJoint *)constraint);
        } else if (constraint->getTypeId() == RigidBodyParticleBallJoint::TYPE_ID) {
            renderRigidBodyParticleBallJoint(*(RigidBodyParticleBallJoint *)constraint);
        } else if (constraint->getTypeId() == RigidBodySpring::TYPE_ID) {
            renderSpring(*(RigidBodySpring *)constraint);
        } else if (constraint->getTypeId() == DistanceJoint::TYPE_ID) {
            renderDistanceJoint(*(DistanceJoint *)constraint);
        } else if (constraint->getTypeId() == DamperJoint::TYPE_ID) {
            renderDamperJoint(*(DamperJoint *)constraint);
        }
    }

    auto *cd = (DistanceFieldCollisionDetection *)Simulation::getCurrent()->getTimeStep()->getCollisionDetection();
    if (cd && (m_renderSDF || m_renderAABB || (m_renderBVHDepth >= 0) || (m_renderBVHDepthTets >= 0))) {
        std::vector<CollisionDetection::CollisionObject *> &collisionObjects = cd->getCollisionObjects();
        for (auto &collisionObject : collisionObjects) {
            if (m_renderAABB) renderAABB(collisionObject->m_aabb);

            if (m_renderSDF) renderSDF(collisionObject);

            if (m_renderBVHDepth >= 0) {
                if (cd->isDistanceFieldCollisionObject(collisionObject)) {
                    const PointCloudBSH &bvh =
                            ((DistanceFieldCollisionDetection::DistanceFieldCollisionObject *)collisionObject)->m_bvh;

                    std::function<bool(unsigned int, unsigned int)> predicate =
                            [&](unsigned int node_index, unsigned int depth) { return (int)depth <= m_renderBVHDepth; };
                    std::function<void(unsigned int, unsigned int)> cb = [&](unsigned int node_index,
                                                                             unsigned int depth) {
                        if (depth == m_renderBVHDepth) {
                            const BoundingSphere &bs = bvh.hull(node_index);
                            if (collisionObject->m_bodyType ==
                                CollisionDetection::CollisionObject::RigidBodyCollisionObjectType) {
                                RigidBody *body = rb[collisionObject->m_bodyIndex];
                                const Vector3r &sphere_x = bs.x();
                                const Vector3r sphere_x_w = body->getRotation() * sphere_x + body->getPosition();
                                MiniGL::drawSphere(sphere_x_w, std::max((float)bs.r(), 0.05f), staticColor);
                            } else
                                MiniGL::drawSphere(bs.x(), std::max((float)bs.r(), 0.05f), staticColor);
                        }
                    };

                    bvh.traverse_depth_first(predicate, cb);
                }
            }

            if (m_renderBVHDepthTets >= 0) {
                if (cd->isDistanceFieldCollisionObject(collisionObject) &&
                    (collisionObject->m_bodyType == CollisionDetection::CollisionObject::TetModelCollisionObjectType)) {
                    TetMeshBSH &bvh = ((DistanceFieldCollisionDetection::DistanceFieldCollisionObject *)collisionObject)
                                              ->m_bvhTets;

                    std::function<bool(unsigned int, unsigned int)> predicate = [&](unsigned int node_index,
                                                                                    unsigned int depth) {
                        return (int)depth <= m_renderBVHDepthTets;
                    };
                    std::function<void(unsigned int, unsigned int)> cb = [&](unsigned int node_index,
                                                                             unsigned int depth) {
                        if (depth == m_renderBVHDepthTets) {
                            const BoundingSphere &bs = bvh.hull(node_index);
                            const Vector3r &sphere_x = bs.x();
                            MiniGL::drawSphere(sphere_x, std::max((float)bs.r(), 0.05f), staticColor);
                        }
                    };

                    bvh.traverse_depth_first(predicate, cb);
                }
            }
        }
    }

    const Vector3r refOffset(0, 0, 0);
    const ParticleData &pd = model->getParticles();
    if (m_renderRefTets || m_renderTets) {
        shaderBegin(surfaceColor);

        for (auto &i : model->getTetModels()) {
            const utility::IndexedTetMesh &mesh = i->getParticleMesh();
            const unsigned int nTets = mesh.numTets();
            const unsigned int *indices = mesh.getTets().data();
            const unsigned int offset = i->getIndexOffset();

            const Vector3r &ix = i->getInitialX();
            const Matrix3r &R = i->getInitialR();

            for (unsigned int j = 0; j < nTets; j++) {
                if (m_renderTets) {
                    const Vector3r &x0 = pd.getPosition(indices[4 * j] + offset);
                    const Vector3r &x1 = pd.getPosition(indices[4 * j + 1] + offset);
                    const Vector3r &x2 = pd.getPosition(indices[4 * j + 2] + offset);
                    const Vector3r &x3 = pd.getPosition(indices[4 * j + 3] + offset);
                    MiniGL::drawTetrahedron(x0, x1, x2, x3, surfaceColor);
                }
                if (m_renderRefTets) {
                    // 					const Vector3r &x0 = R.transpose() * (pd.getPosition0(indices[4 * j + 0]
                    // + offset)
                    // - ix); 					const Vector3r &x1 = R.transpose() *
                    // (pd.getPosition0(indices[4 * j + 1] + offset) - ix); 					const Vector3r &x2 =
                    // R.transpose()
                    // * (pd.getPosition0(indices[4 * j + 2] + offset) - ix); 					const
                    // Vector3r &x3 = R.transpose() * (pd.getPosition0(indices[4 * j + 3] + offset) - ix);
                    const Vector3r &x0 = pd.getPosition0(indices[4 * j] + offset) + refOffset;
                    const Vector3r &x1 = pd.getPosition0(indices[4 * j + 1] + offset) + refOffset;
                    const Vector3r &x2 = pd.getPosition0(indices[4 * j + 2] + offset) + refOffset;
                    const Vector3r &x3 = pd.getPosition0(indices[4 * j + 3] + offset) + refOffset;

                    MiniGL::drawTetrahedron(x0, x1, x2, x3, staticColor);
                }
            }
        }
        shaderEnd();
    }

    float red[4] = {0.8f, 0.0f, 0.0f, 1};
    for (unsigned int m_selectedParticle : m_selectedParticles) {
        MiniGL::drawSphere(pd.getPosition(m_selectedParticle), 0.08f, red);
    }

    MiniGL::drawTime(TimeManager::getCurrent()->getTime());
}

void DemoBase::renderTriangleModels() {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const ParticleData &pd = model->getParticles();
    float surfaceColor[4] = {0.8f, 0.9f, 0.2f, 1};

    shaderTexBegin(surfaceColor);

    for (auto &i : model->getTriangleModels()) {
        // mesh
        const utility::IndexedFaceMesh &mesh = i->getParticleMesh();
        const unsigned int offset = i->getIndexOffset();
        Visualization::drawTexturedMesh(pd, mesh, offset, surfaceColor);
    }

    shaderTexEnd();
}

void DemoBase::renderTetModels() {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const ParticleData &pd = model->getParticles();
    float surfaceColor[4] = {0.1f, 0.4f, 0.7f, 1};

    shaderBegin(surfaceColor);

    for (auto &i : model->getTetModels()) {
        const VertexData &vdVis = i->getVisVertices();
        if (vdVis.size() > 0) {
            const utility::IndexedFaceMesh &visMesh = i->getVisMesh();
            Visualization::drawMesh(vdVis, visMesh, 0, surfaceColor);
        } else {
            const utility::IndexedFaceMesh &surfaceMesh = i->getSurfaceMesh();
            const unsigned int offset = i->getIndexOffset();
            Visualization::drawMesh(pd, surfaceMesh, offset, surfaceColor);
        }
    }

    shaderEnd();
}

void DemoBase::renderAABB(AABB &aabb) {
    Vector3r p1, p2;
    glBegin(GL_LINES);
    for (unsigned char i = 0; i < 12; i++) {
        AABB::getEdge(aabb, i, p1, p2);
        glVertex3d(p1[0], p1[1], p1[2]);
        glVertex3d(p2[0], p2[1], p2[2]);
    }
    glEnd();
}

void DemoBase::renderSDF(CollisionDetection::CollisionObject *co) {
    auto *cd = (DistanceFieldCollisionDetection *)Simulation::getCurrent()->getTimeStep()->getCollisionDetection();
    if ((!cd->isDistanceFieldCollisionObject(co)) ||
        (co->m_bodyType != CollisionDetection::CollisionObject::RigidBodyCollisionObjectType))
        return;

    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[co->m_bodyIndex];

    const Vector3r &com = rb->getPosition();
    const Matrix3r &R = rb->getTransformationR();
    const Vector3r &v1 = rb->getTransformationV1();
    const Vector3r &v2 = rb->getTransformationV2();

    auto *dfco = (DistanceFieldCollisionDetection::DistanceFieldCollisionObject *)co;
    const Vector3r &startX = co->m_aabb.m_p[0];
    const Vector3r &endX = co->m_aabb.m_p[1];
    Vector3r diff = endX - startX;
    const unsigned int steps = 20;
    Vector3r stepSize = (1.0 / steps) * diff;
    for (Real x = startX[0]; x < endX[0]; x += stepSize[0]) {
        for (Real y = startX[1]; y < endX[1]; y += stepSize[1]) {
            for (Real z = startX[2]; z < endX[2]; z += stepSize[2]) {
                Vector3r pos_w(x, y, z);
                const Vector3r pos = R * (pos_w - com) + v1;
                const double dist = dfco->distance(pos.template cast<double>(), 0.0);

                if (dist < 0.0) {
                    float col[4] = {(float)-dist, 0.0f, 0.0f, 1.0f};
                    MiniGL::drawPoint(pos_w, 3.0f, col);
                }
            }
        }
    }
}

void DemoBase::renderBallJoint(BallJoint &bj) { MiniGL::drawSphere(bj.m_jointInfo.col(2), 0.15f, m_jointColor); }

void DemoBase::renderRigidBodyParticleBallJoint(RigidBodyParticleBallJoint &bj) {
    MiniGL::drawSphere(bj.m_jointInfo.col(1), 0.1f, m_jointColor);
}

void DemoBase::renderBallOnLineJoint(BallOnLineJoint &bj) {
    MiniGL::drawSphere(bj.m_jointInfo.col(5), 0.1f, m_jointColor);
    MiniGL::drawCylinder(bj.m_jointInfo.col(5) - bj.m_jointInfo.col(7), bj.m_jointInfo.col(5) + bj.m_jointInfo.col(7),
                         m_jointColor, 0.05f);
}

void DemoBase::renderHingeJoint(HingeJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 4);
    const Vector3r &axis_local = joint.m_jointInfo.block<3, 1>(0, 6);
    const Vector3r axis = rb->getRotation().matrix() * axis_local;

    MiniGL::drawSphere(c - 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawSphere(c + 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - 0.5 * axis, c + 0.5 * axis, m_jointColor, 0.05f);
}

void DemoBase::renderUniversalJoint(UniversalJoint &uj) {
    MiniGL::drawSphere(uj.m_jointInfo.col(4) - 0.5 * uj.m_jointInfo.col(6), 0.1f, m_jointColor);
    MiniGL::drawSphere(uj.m_jointInfo.col(4) + 0.5 * uj.m_jointInfo.col(6), 0.1f, m_jointColor);
    MiniGL::drawSphere(uj.m_jointInfo.col(5) - 0.5 * uj.m_jointInfo.col(7), 0.1f, m_jointColor);
    MiniGL::drawSphere(uj.m_jointInfo.col(5) + 0.5 * uj.m_jointInfo.col(7), 0.1f, m_jointColor);
    MiniGL::drawCylinder(uj.m_jointInfo.col(4) - 0.5 * uj.m_jointInfo.col(6),
                         uj.m_jointInfo.col(4) + 0.5 * uj.m_jointInfo.col(6), m_jointColor, 0.05f);
    MiniGL::drawCylinder(uj.m_jointInfo.col(5) - 0.5 * uj.m_jointInfo.col(7),
                         uj.m_jointInfo.col(5) + 0.5 * uj.m_jointInfo.col(7), m_jointColor, 0.05f);
}

void DemoBase::renderSliderJoint(SliderJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    Quaternionr qR0;
    qR0.coeffs() = joint.m_jointInfo.col(1);
    const Vector3r &c = rb->getPosition();
    Vector3r axis = qR0.matrix().col(0);
    MiniGL::drawSphere(c, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void DemoBase::renderTargetPositionMotorSliderJoint(TargetPositionMotorSliderJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    const Vector3r &c = rb->getPosition();
    Vector3r axis = joint.m_jointInfo.block<3, 1>(0, 1);
    MiniGL::drawSphere(c, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void DemoBase::renderTargetVelocityMotorSliderJoint(TargetVelocityMotorSliderJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    Quaternionr qR0;
    qR0.coeffs() = joint.m_jointInfo.col(1);
    const Vector3r &c = rb->getPosition();
    Vector3r axis = qR0.matrix().col(0);
    MiniGL::drawSphere(c, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void DemoBase::renderTargetAngleMotorHingeJoint(TargetAngleMotorHingeJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 5);
    const Vector3r &axis_local = joint.m_jointInfo.block<3, 1>(0, 7);
    const Vector3r axis = rb->getRotation().matrix() * axis_local;

    MiniGL::drawSphere(c - 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawSphere(c + 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - 0.5 * axis, c + 0.5 * axis, m_jointColor, 0.05f);
}

void DemoBase::renderTargetVelocityMotorHingeJoint(TargetVelocityMotorHingeJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();

    const Vector3r &c = joint.m_jointInfo.block<3, 1>(0, 5);
    const Vector3r axis = joint.m_jointInfo.block<3, 1>(0, 7);

    MiniGL::drawSphere(c - 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawSphere(c + 0.5 * axis, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - 0.5 * axis, c + 0.5 * axis, m_jointColor, 0.05f);
}

void DemoBase::renderRigidBodyContact(RigidBodyContactConstraint &cc) {
    float col1[4] = {0.0f, 0.6f, 0.2f, 1};
    float col2[4] = {0.6f, 0.0f, 0.2f, 1};
    MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
    MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
    MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void DemoBase::renderParticleRigidBodyContact(ParticleRigidBodyContactConstraint &cc) {
    float col1[4] = {0.0f, 0.6f, 0.2f, 1};
    float col2[4] = {0.6f, 0.0f, 0.2f, 1};
    MiniGL::drawPoint(cc.m_constraintInfo.col(0), 5.0f, col1);
    MiniGL::drawPoint(cc.m_constraintInfo.col(1), 5.0f, col2);
    MiniGL::drawVector(cc.m_constraintInfo.col(1), cc.m_constraintInfo.col(1) + cc.m_constraintInfo.col(2), 1.0f, col2);
}

void DemoBase::renderSpring(RigidBodySpring &s) {
    MiniGL::drawSphere(s.m_jointInfo.col(2), 0.1f, m_jointColor);
    MiniGL::drawSphere(s.m_jointInfo.col(3), 0.1f, m_jointColor);
    MiniGL::drawCylinder(s.m_jointInfo.col(2), s.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void DemoBase::renderDistanceJoint(DistanceJoint &j) {
    MiniGL::drawSphere(j.m_jointInfo.col(2), 0.1f, m_jointColor);
    MiniGL::drawSphere(j.m_jointInfo.col(3), 0.1f, m_jointColor);
    MiniGL::drawCylinder(j.m_jointInfo.col(2), j.m_jointInfo.col(3), m_jointColor, 0.05f);
}

void DemoBase::renderDamperJoint(DamperJoint &joint) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const SimulationModel::RigidBodyVector &rigidBodies = model->getRigidBodies();
    RigidBody *rb = rigidBodies[joint.m_bodies[0]];

    Quaternionr qR0;
    qR0.coeffs() = joint.m_jointInfo.col(1);
    const Vector3r &c = rb->getPosition();
    Vector3r axis = qR0.matrix().col(0);
    MiniGL::drawSphere(c, 0.1f, m_jointColor);
    MiniGL::drawCylinder(c - axis, c + axis, m_jointColor, 0.05f);
}

void DemoBase::mouseMove(int x, int y, void *clientData) {
    auto *base = (DemoBase *)clientData;
    SimulationModel *model = Simulation::getCurrent()->getModel();

    Vector3r mousePos;
    MiniGL::unproject(x, y, mousePos);
    const Vector3r diff = mousePos - base->m_oldMousePos;

    TimeManager *tm = TimeManager::getCurrent();
    const Real h = tm->getTimeStepSize();

    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    for (unsigned int m_selectedBodie : base->m_selectedBodies) {
        const Real mass = rb[m_selectedBodie]->getMass();
        if (mass != 0.0) rb[m_selectedBodie]->getVelocity() += 3.0 / h * diff;
    }
    ParticleData &pd = model->getParticles();
    for (unsigned int m_selectedParticle : base->m_selectedParticles) {
        const Real mass = pd.getMass(m_selectedParticle);
        if (mass != 0.0) pd.getVelocity(m_selectedParticle) += 5.0 * diff / h;
    }
    base->m_oldMousePos = mousePos;
}

void DemoBase::selection(const Vector2i &start, const Vector2i &end, void *clientData) {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    auto *base = (DemoBase *)clientData;

    std::vector<unsigned int> hits;

    base->m_selectedParticles.clear();
    ParticleData &pd = model->getParticles();
    if (pd.size() > 0)
        Selection::selectRect(start, end, &pd.getPosition(0), &pd.getPosition(pd.size() - 1),
                              base->m_selectedParticles);

    base->m_selectedBodies.clear();
    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> x;
    x.resize(rb.size());
    for (unsigned int i = 0; i < rb.size(); i++) {
        x[i] = rb[i]->getPosition();
    }

    if (!rb.empty()) Selection::selectRect(start, end, &x[0], &x[rb.size() - 1], base->m_selectedBodies);
    if ((!base->m_selectedBodies.empty()) || (!base->m_selectedParticles.empty()))
        MiniGL::setMouseMoveFunc(2, mouseMove);
    else
        MiniGL::setMouseMoveFunc(-1, nullptr);

    MiniGL::unproject(end[0], end[1], base->m_oldMousePos);
}

void DemoBase::reset() {}