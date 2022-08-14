//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/gui/opengl/mini_gl.h"
#include "vox.sph/gui/opengl/shader.h"
#include "vox.sph/gui/tweakbar/simulator_gui_tweakbar.h"
#include "vox.sph/position_based_dynamics_wrapper/pbd_wrapper.h"

namespace vox {
class PBD_Simulator_GUI_TweakBar : public Simulator_GUI_TweakBar {
protected:
    PBDWrapper *m_pbdWrapper;
    bool m_drawAABB;
    int m_drawBVHDepth;
    bool m_drawSDF;
    Shader *m_shader;
    float m_jointColor[4];

    Shader *createShader(const std::string &vertexShader,
                         const std::string &geometryShader,
                         const std::string &fragmentShader);
    void initShader();

    void shaderBegin(const float *col);
    void shaderEnd();

    void renderAABB(AABB &aabb);
    void renderBVH();
    void renderSDF();
    void renderSDF(CollisionObject *co);
    void renderTriangleModels();
    void renderTetModels();
    void renderConstraints();
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
    void renderDamperJoint(DamperJoint &joint);

    static void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData);
    static void TW_CALL getVelocityUpdateMethod(void *value, void *clientData);
    static void TW_CALL setDampingCoeff(const void *value, void *clientData);
    static void TW_CALL getDampingCoeff(void *value, void *clientData);
    static void TW_CALL setMaxIterations(const void *value, void *clientData);
    static void TW_CALL getMaxIterations(void *value, void *clientData);
    static void TW_CALL setMaxIterationsV(const void *value, void *clientData);
    static void TW_CALL getMaxIterationsV(void *value, void *clientData);
    static void TW_CALL setContactTolerance(const void *value, void *clientData);
    static void TW_CALL getContactTolerance(void *value, void *clientData);
    static void TW_CALL setContactStiffnessRigidBody(const void *value, void *clientData);
    static void TW_CALL getContactStiffnessRigidBody(void *value, void *clientData);
    static void TW_CALL setContactStiffnessParticleRigidBody(const void *value, void *clientData);
    static void TW_CALL getContactStiffnessParticleRigidBody(void *value, void *clientData);
    static void TW_CALL setStiffness(const void *value, void *clientData);
    static void TW_CALL getStiffness(void *value, void *clientData);
    static void TW_CALL setXXStiffness(const void *value, void *clientData);
    static void TW_CALL getXXStiffness(void *value, void *clientData);
    static void TW_CALL setYYStiffness(const void *value, void *clientData);
    static void TW_CALL getYYStiffness(void *value, void *clientData);
    static void TW_CALL setXYStiffness(const void *value, void *clientData);
    static void TW_CALL getXYStiffness(void *value, void *clientData);
    static void TW_CALL setXYPoissonRatio(const void *value, void *clientData);
    static void TW_CALL getXYPoissonRatio(void *value, void *clientData);
    static void TW_CALL setYXPoissonRatio(const void *value, void *clientData);
    static void TW_CALL getYXPoissonRatio(void *value, void *clientData);
    static void TW_CALL setNormalizeStretch(const void *value, void *clientData);
    static void TW_CALL getNormalizeStretch(void *value, void *clientData);
    static void TW_CALL setNormalizeShear(const void *value, void *clientData);
    static void TW_CALL getNormalizeShear(void *value, void *clientData);
    static void TW_CALL setBendingStiffness(const void *value, void *clientData);
    static void TW_CALL getBendingStiffness(void *value, void *clientData);
    static void TW_CALL setBendingMethod(const void *value, void *clientData);
    static void TW_CALL getBendingMethod(void *value, void *clientData);
    static void TW_CALL setClothSimulationMethod(const void *value, void *clientData);
    static void TW_CALL getClothSimulationMethod(void *value, void *clientData);
    static void TW_CALL setSolidSimulationMethod(const void *value, void *clientData);
    static void TW_CALL getSolidSimulationMethod(void *value, void *clientData);

public:
    PBD_Simulator_GUI_TweakBar(SimulatorBase *base, PBDWrapper *pbdWrapper);
    virtual ~PBD_Simulator_GUI_TweakBar();

public:
    virtual void init(int argc, char **argv, const char *name);
    virtual void render();
    virtual void initSimulationParameterGUI();

public:
    template <class PositionData>
    static void drawMesh(const PositionData &pd,
                         const utility::IndexedFaceMesh &mesh,
                         const unsigned int offset,
                         const float *const color);
};

template <class PositionData>
void PBD_Simulator_GUI_TweakBar::drawMesh(const PositionData &pd,
                                          const utility::IndexedFaceMesh &mesh,
                                          const unsigned int offset,
                                          const float *const color) {
    // draw mesh
    const unsigned int *faces = mesh.getFaces().data();
    const unsigned int nFaces = mesh.numFaces();
    const Vector3r *vertexNormals = mesh.getVertexNormals().data();

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &pd.getPosition(offset)[0]);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_REAL, GL_FALSE, 0, &vertexNormals[0][0]);
    } else {
        float speccolor[4] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f);
        glColor3fv(color);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_REAL, 0, &pd.getPosition(0)[0]);
        glNormalPointer(GL_REAL, 0, &vertexNormals[0][0]);
    }

    glDrawElements(GL_TRIANGLES, (GLsizei)3 * mesh.numFaces(), GL_UNSIGNED_INT, mesh.getFaces().data());

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(2);
    } else {
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
}

}  // namespace vox