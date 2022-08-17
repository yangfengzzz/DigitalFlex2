//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <cmath>
#include <iostream>

#include "apps.pbd/common/demo_base.h"
#include "apps.pbd/common/mini_gl.h"
#include "apps.pbd/common/tweakbar_parameters.h"
#include "vox.base/file_system.h"
#include "vox.base/logging.h"
#include "vox.base/obj_loader.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.pbd/simulation.h"
#include "vox.pbd/time_step_controller.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace vox;
using namespace Eigen;
using namespace std;

void initParameters();
void timeStep();
void buildModel();
void createBodyModel();
void render();
void reset();

DemoBase *base;

const Real width = static_cast<Real>(0.4);
const Real height = static_cast<Real>(0.4);
const Real depth = static_cast<Real>(2.0);

// main
int main(int argc, char **argv) {
    REPORT_MEMORY_LEAKS

    base = new DemoBase();
    base->init(argc, argv, "Rigid body demo");

    auto *model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    initParameters();

    Simulation::getCurrent()->getTimeStep()->setValue<unsigned int>(TimeStepController::NUM_SUB_STEPS, 1);
    Simulation::getCurrent()->getTimeStep()->setValue<unsigned int>(TimeStepController::MAX_ITERATIONS, 5);

    // OpenGL
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(60.0, 0.1f, 500.0, Vector3r(6.0, 0.0, 18.0), Vector3r(6.0, -4.0, 0.0));

    MiniGL::mainLoop();

    Timing::printAverageTimes();
    Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

void initParameters() {
    TwRemoveAllVars(MiniGL::getTweakBar());
    TweakBarParameters::cleanup();

    MiniGL::initTweakBarParameters();

    TweakBarParameters::createParameterGUI();
    TweakBarParameters::createParameterObjectGUI(base);
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getModel());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getTimeStep());
}

void reset() {
    Timing::printAverageTimes();
    Timing::reset();

    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();
}

void timeStep() {
    const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
    if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime())) base->setValue(DemoBase::PAUSE, true);

    if (base->getValue<bool>(DemoBase::PAUSE)) return;

    // Simulation code
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const auto numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
    for (unsigned int i = 0; i < numSteps; i++) {
        // set target angle of motors for an animation
        const Real currentTargetAngle =
                static_cast<Real>(M_PI * 0.5) -
                static_cast<Real>(M_PI * 0.25) *
                        (cos(static_cast<Real>(0.25) * TimeManager::getCurrent()->getTime()) + 1.0);
        SimulationModel::ConstraintVector &constraints = model->getConstraints();
        TargetAngleMotorHingeJoint &joint1 = (*(TargetAngleMotorHingeJoint *)constraints[8]);
        TargetVelocityMotorHingeJoint &joint2 = (*(TargetVelocityMotorHingeJoint *)constraints[10]);
        joint1.setTarget(currentTargetAngle);
        joint2.setTarget(3.5);

        const Real currentTargetPos =
                static_cast<Real>(1.5) * sin(static_cast<Real>(2.0) * TimeManager::getCurrent()->getTime());
        TargetPositionMotorSliderJoint &joint3 = (*(TargetPositionMotorSliderJoint *)constraints[12]);
        joint3.setTarget(currentTargetPos);

        Real currentTargetVel = 0.25;
        if (((int)(0.25 * TimeManager::getCurrent()->getTime())) % 2 == 1) currentTargetVel = -currentTargetVel;
        TargetVelocityMotorSliderJoint &joint4 = (*(TargetVelocityMotorSliderJoint *)constraints[14]);
        joint4.setTarget(currentTargetVel);

        START_TIMING("SimStep")
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG
    }
}

void buildModel() {
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));
    createBodyModel();

    std::cout << "Joint types in simulation: \n";
    std::cout << "--------------------------\n\n";
    std::cout << "row 1:    ball joint,                  ball-on-line joint,             hinge joint,               "
                 "universal joint\n";
    std::cout << "row 2:    target angle hinge motor,    target velocity hinge motor,    target position motor,     "
                 "target velocity motor\n";
    std::cout << "row 3:    spring,                      distance joint,                 slider joint\n";
}

void render() {
    base->render();

    // float textColor[4] = { 0.0, .2f, .4f, 1 };
    // MiniGL::drawStrokeText(-0.5, 1.5, 1.0, 0.002f, "ball joint", 11, textColor);
    // MiniGL::drawStrokeText(3.0, 1.5, 1.0, 0.002f, "ball-on-line joint", 19, textColor);
    // MiniGL::drawStrokeText(7.3f, 1.5, 1.0, 0.002f, "hinge joint", 12, textColor);
    // MiniGL::drawStrokeText(11.2f, 1.5, 1.0, 0.002f, "universal joint", 15, textColor);

    // MiniGL::drawStrokeText(-1.0, -4.0, 1.0, 0.002f, "motor hinge joint", 17, textColor);
    // MiniGL::drawStrokeText(3.4f, -4.0, 1.0, 0.002f, "slider joint", 12, textColor);
    // MiniGL::drawStrokeText(6.6f, -4.0, 1.0, 0.002f, "target position motor", 21, textColor);
    // MiniGL::drawStrokeText(10.6f, -4.0, 1.0, 0.002f, "target velocity motor", 21, textColor);

    // MiniGL::drawStrokeText(-0.25, -9.5, 1.0, 0.002f, "spring", 6, textColor);
    // MiniGL::drawStrokeText(3.3, -9.5, 1.0, 0.002f, "distance joint", 14, textColor);

    MiniGL::drawTime(TimeManager::getCurrent()->getTime());
}

// Compute diagonal inertia tensor
Vector3r computeInertiaTensorBox(const Real mass, const Real width, const Real height, const Real depth) {
    const Real Ix = (mass / static_cast<Real>(12.0)) * (height * height + depth * depth);
    const Real Iy = (mass / static_cast<Real>(12.0)) * (width * width + depth * depth);
    const Real Iz = (mass / static_cast<Real>(12.0)) * (width * width + height * height);
    return Vector3r(Ix, Iy, Iz);
}

void loadObj(const std::string &filename, VertexData &vd, utility::IndexedFaceMesh &mesh, const Vector3r &scale) {
    std::vector<utility::OBJLoader::Vec3f> x;
    std::vector<utility::OBJLoader::Vec3f> normals;
    std::vector<utility::OBJLoader::Vec2f> texCoords;
    std::vector<utility::MeshFaceIndices> faces;
    utility::OBJLoader::Vec3f s = {(float)scale[0], (float)scale[1], (float)scale[2]};
    utility::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

    mesh.release();
    const auto nPoints = (unsigned int)x.size();
    const auto nFaces = (unsigned int)faces.size();
    const auto nTexCoords = (unsigned int)texCoords.size();
    mesh.initMesh(nPoints, nFaces * 2, nFaces);
    vd.reserve(nPoints);
    for (unsigned int i = 0; i < nPoints; i++) {
        vd.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
    }
    for (unsigned int i = 0; i < nTexCoords; i++) {
        mesh.addUV(texCoords[i][0], texCoords[i][1]);
    }
    for (unsigned int i = 0; i < nFaces; i++) {
        // Reduce the indices by one
        int posIndices[3];
        int texIndices[3];
        for (int j = 0; j < 3; j++) {
            posIndices[j] = faces[i].posIndices[j] - 1;
            if (nTexCoords > 0) {
                texIndices[j] = faces[i].texIndices[j] - 1;
                mesh.addUVIndex(texIndices[j]);
            }
        }

        mesh.addFace(&posIndices[0]);
    }
    mesh.buildNeighbors();

    mesh.updateNormals(vd, 0);
    mesh.updateVertexNormals(vd);

    LOGI("Number of triangles: {}", nFaces)
    LOGI("Number of vertices: {}", nPoints)
}

/** Create the rigid body model
 */
void createBodyModel() {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();

    string fileName = utility::FileSystem::normalizePath(base->getExePath() + "/resources/models/cube.obj");
    utility::IndexedFaceMesh mesh;
    VertexData vd;
    loadObj(fileName, vd, mesh, Vector3r(width, height, depth));
    mesh.setFlatShading(true);
    utility::IndexedFaceMesh meshStatic;
    VertexData vdStatic;
    loadObj(fileName, vdStatic, meshStatic, Vector3r(0.5, 0.5, 0.5));
    meshStatic.setFlatShading(true);

    // static body
    const unsigned int numberOfBodies = 33;
    rb.resize(numberOfBodies);
    Real startX = 0.0;
    Real startY = 6.5;
    for (unsigned int i = 0; i < 11; i++) {
        if (i % 4 == 0) {
            startY -= 5.5;
            startX = 0.0;
        }

        rb[3 * i] = new RigidBody();
        rb[3 * i]->initBody(0.0, Vector3r(startX, startY, 1.0), computeInertiaTensorBox(1.0, 0.5, 0.5, 0.5),
                            Quaternionr(1.0, 0.0, 0.0, 0.0), vdStatic, meshStatic);

        // dynamic body
        rb[3 * i + 1] = new RigidBody();
        rb[3 * i + 1]->initBody(1.0, Vector3r(startX, startY - static_cast<Real>(0.25), static_cast<Real>(2.0)),
                                computeInertiaTensorBox(1.0, width, height, depth), Quaternionr(1.0, 0.0, 0.0, 0.0), vd,
                                mesh);

        // dynamic body
        rb[3 * i + 2] = new RigidBody();
        rb[3 * i + 2]->initBody(1.0, Vector3r(startX, startY - static_cast<Real>(0.25), static_cast<Real>(4.0)),
                                computeInertiaTensorBox(1.0, width, height, depth), Quaternionr(1.0, 0.0, 0.0, 0.0), vd,
                                mesh);

        startX += 4.0;
    }

    Real jointY = 0.75;
    model->addBallJoint(0, 1, Vector3r(0.25, jointY, 1.0));
    model->addBallJoint(1, 2, Vector3r(0.25, jointY, 3.0));

    model->addBallOnLineJoint(3, 4, Vector3r(4.25, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(4, 5, Vector3r(4.25, jointY, 3.0));

    model->addHingeJoint(6, 7, Vector3r(8.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(7, 8, Vector3r(8.25, jointY, 3.0));

    model->addUniversalJoint(9, 10, Vector3r(12.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0), Vector3r(0.0, 1.0, 0.0));
    model->addBallJoint(10, 11, Vector3r(12.25, jointY, 3.0));

    jointY -= 5.5;
    model->addTargetAngleMotorHingeJoint(12, 13, Vector3r(0.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(13, 14, Vector3r(0.25, jointY, 3.0));

    model->addTargetVelocityMotorHingeJoint(15, 16, Vector3r(4.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(16, 17, Vector3r(4.25, jointY, 3.0));

    model->addTargetPositionMotorSliderJoint(18, 19, Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(19, 20, Vector3r(8.25, jointY, 3.0));

    model->addTargetVelocityMotorSliderJoint(21, 22, Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(22, 23, Vector3r(12.25, jointY, 3.0));

    jointY -= 5.5;
    model->addRigidBodySpring(24, 25, Vector3r(0.25, jointY, 1.0), Vector3r(0.25, jointY, 1.0), 50.0);
    model->addBallJoint(25, 26, Vector3r(0.25, jointY, 3.0));

    model->addDistanceJoint(27, 28, Vector3r(4.25, jointY, 1.0), Vector3r(4.25, jointY, 2.0));
    model->addBallJoint(28, 29, Vector3r(4.25, jointY, 3.0));

    model->addSliderJoint(30, 31, Vector3r(1.0, 0.0, 0.0));
    model->addBallJoint(31, 32, Vector3r(8.25, jointY, 3.0));
}
