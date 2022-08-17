//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <iostream>

#include "apps.pbd/common/demo_base.h"
#include "apps.pbd/common/mini_gl.h"
#include "apps.pbd/common/tweakbar_parameters.h"
#include "apps.pbd/generic_constraints/generic_constraints_model.h"
#include "vox.base/file_system.h"
#include "vox.base/logging.h"
#include "vox.base/obj_loader.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.pbd/simulation.h"

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
    base->init(argc, argv, "Generic rigid body demo");

    auto *model = new GenericConstraintsModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    initParameters();

    // OpenGL
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(60.0, 0.1f, 500.0, Vector3r(2.0, 4.0, 15.0), Vector3r(2.0, -2.0, 0.0));

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
        START_TIMING("SimStep")
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG
    }
}

void buildModel() {
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

    createBodyModel();
}

void render() { base->render(); }

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

void createBodyModel() {
    auto *model = (GenericConstraintsModel *)Simulation::getCurrent()->getModel();
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
    const unsigned int numberOfBodies = 6;
    rb.resize(numberOfBodies);
    Real startX = 0.0;
    Real startY = 1.0;
    for (unsigned int i = 0; i < 2; i++) {
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

        if (i == 3) {
            startY -= 5.5;
            startX = 0.0;
        }
    }

    Real jointY = 0.75;
    model->addGenericHingeJoint(0, 1, Vector3r(0.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addGenericBallJoint(1, 2, Vector3r(0.25, jointY, 3.0));

    model->addGenericSliderJoint(3, 4, Vector3r(4.0, jointY, 1.0), Vector3r(1.0, 0.0, 0.0));
    model->addGenericBallJoint(4, 5, Vector3r(4.25, jointY, 3.0));
}
