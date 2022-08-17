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
void createBunnyRodModel();
void render();
void reset();

DemoBase *base;

// bunny rod scene
const int numberOfBodies = 10;
const Real width = static_cast<Real>(1.0);
const Real height = static_cast<Real>(0.1);
const Real depth = static_cast<Real>(0.1);
const Real youngsModulus = static_cast<Real>(209e9);
const Real torsionModulus = static_cast<Real>(79e9);
const Real density = 7800.;

const Real bunnyDensity = 500.;

// main
int main(int argc, char **argv) {
    REPORT_MEMORY_LEAKS

    base = new DemoBase();
    base->init(argc, argv, "Stretch-bending-twisting demo");

    auto *model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    initParameters();

    // OpenGL
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(40.0, 0.1f, 500.0, Vector3r(5.0, 10.0, 30.0), Vector3r(5.0, 0.0, 0.0));

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

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
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
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.01));

    createBunnyRodModel();
}

void render() { base->render(); }

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

/** Create the bunny rod body model
 */
void createBunnyRodModel() {
    SimulationModel *model = Simulation::getCurrent()->getModel();
    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    SimulationModel::ConstraintVector &constraints = model->getConstraints();

    string fileName = utility::FileSystem::normalizePath(base->getExePath() + "/resources/models/cube.obj");
    utility::IndexedFaceMesh mesh;
    VertexData vd;
    loadObj(fileName, vd, mesh, Vector3r(height, width, depth));
    mesh.setFlatShading(true);

    string fileName2 = utility::FileSystem::normalizePath(base->getExePath() + "/resources/models/bunny_10k.obj");
    utility::IndexedFaceMesh mesh2;
    VertexData vd2;
    loadObj(fileName2, vd2, mesh2, Vector3r(2.0, 2.0, 2.0));

    rb.resize(numberOfBodies);
    for (unsigned int i = 0; i < numberOfBodies - 1; i++) {
        rb[i] = new RigidBody();

        Real mass(static_cast<Real>(0.25 * M_PI) * width * height * depth * density);

        const Real Iy = static_cast<Real>(1. / 12.) * mass *
                        (static_cast<Real>(3.) * (static_cast<Real>(0.25) * height * height) + width * width);
        const Real Iz = static_cast<Real>(1. / 12.) * mass *
                        (static_cast<Real>(3.) * (static_cast<Real>(0.25) * depth * depth) + width * width);
        const Real Ix = static_cast<Real>(0.25) * mass * (static_cast<Real>(0.25) * (height * height + depth * depth));
        Vector3r inertia(Iy, Ix, Iz);  // rod axis along y-axis

        Real angle(static_cast<Real>(M_PI_2));
        Vector3r axis(0., 0., 1.);
        AngleAxisr angleAxis(angle, axis);
        Quaternionr rotation(angleAxis);

        rb[i]->initBody(mass, Vector3r((Real)i * width, 0.0, 0.0), inertia, rotation, vd, mesh);
    }
    // Make first body static
    rb[0]->setMass(0.0);

    // bunny
    const Quaternionr q(AngleAxisr(static_cast<Real>(1.0 / 6.0 * M_PI), Vector3r(0.0, 0.0, 1.0)));
    const Vector3r t(static_cast<Real>(0.411) + (static_cast<Real>(numberOfBodies - 1.0)) * width,
                     -static_cast<Real>(1.776), static_cast<Real>(0.356));
    rb[numberOfBodies - 1] = new RigidBody();
    rb[numberOfBodies - 1]->initBody(bunnyDensity, t, q, vd2, mesh2);

    constraints.reserve(numberOfBodies - 1);
    for (unsigned int i = 0; i < numberOfBodies - 2; i++) {
        model->addStretchBendingTwistingConstraint(
                i, i + 1, Vector3r(static_cast<Real>(i) * width + static_cast<Real>(0.5) * width, 0.0, 0.0),
                0.25 * (height + depth), width, youngsModulus, torsionModulus);
    }
    unsigned int i(numberOfBodies - 2);
    model->addBallJoint(i, i + 1, Vector3r(static_cast<Real>(i) * width + static_cast<Real>(0.5) * width, 0.0, 0.0));

    Simulation::getCurrent()->getTimeStep()->setValue(TimeStepController::MAX_ITERATIONS, 200);
}
