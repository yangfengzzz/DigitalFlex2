//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <Eigen/Dense>
#include <iostream>

#include "apps.pbd/common/demo_base.h"
#include "apps.pbd/common/mini_gl.h"
#include "apps.pbd/common/tweakbar_parameters.h"
#include "vox.base/file_system.h"
#include "vox.base/logging.h"
#include "vox.base/obj_loader.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.pbd/distance_field_collision_detection.h"
#include "vox.pbd/simulation.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace vox;
using namespace Eigen;
using namespace std;

void timeStep();
void buildModel();
void createMesh();
void render();
void reset();
void TW_CALL setBendingMethod(const void* value, void* clientData);
void TW_CALL getBendingMethod(void* value, void* clientData);
void TW_CALL setSimulationMethod(const void* value, void* clientData);
void TW_CALL getSimulationMethod(void* value, void* clientData);
void TW_CALL setBendingStiffness(const void* value, void* clientData);
void TW_CALL getBendingStiffness(void* value, void* clientData);
void TW_CALL setDistanceStiffness(const void* value, void* clientData);
void TW_CALL getDistanceStiffness(void* value, void* clientData);
void TW_CALL setXXStiffness(const void* value, void* clientData);
void TW_CALL getXXStiffness(void* value, void* clientData);
void TW_CALL setYYStiffness(const void* value, void* clientData);
void TW_CALL getYYStiffness(void* value, void* clientData);
void TW_CALL setXYStiffness(const void* value, void* clientData);
void TW_CALL getXYStiffness(void* value, void* clientData);
void TW_CALL setXYPoissonRatio(const void* value, void* clientData);
void TW_CALL getXYPoissonRatio(void* value, void* clientData);
void TW_CALL setYXPoissonRatio(const void* value, void* clientData);
void TW_CALL getYXPoissonRatio(void* value, void* clientData);
void TW_CALL setNormalizeStretch(const void* value, void* clientData);
void TW_CALL getNormalizeStretch(void* value, void* clientData);
void TW_CALL setNormalizeShear(const void* value, void* clientData);
void TW_CALL getNormalizeShear(void* value, void* clientData);

const int nRows = 50;
const int nCols = 50;
const Real width = 10.0;
const Real height = 10.0;
short simulationMethod = 2;
short bendingMethod = 2;
Real distanceStiffness = 1.0;
Real xxStiffness = 1.0;
Real yyStiffness = 1.0;
Real xyStiffness = 1.0;
Real xyPoissonRatio = 0.3;
Real yxPoissonRatio = 0.3;
bool normalizeStretch = false;
bool normalizeShear = false;
Real bendingStiffness = 0.01;
bool doPause = true;
DemoBase* base;
DistanceFieldCollisionDetection cd;

// main
int main(int argc, char** argv) {
    REPORT_MEMORY_LEAKS

    base = new DemoBase();
    base->init(argc, argv, "Cloth demo");

    auto* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    base->createParameterGUI();

    // OpenGL
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3r(0.0, 10.0, 25.0), Vector3r(0.0, 0.0, 0.0));

    TwType enumType2 = TwDefineEnum("SimulationMethodType", nullptr, 0);
    TwAddVarCB(MiniGL::getTweakBar(), "SimulationMethod", enumType2, setSimulationMethod, getSimulationMethod,
               &simulationMethod,
               " label='Simulation method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain "
               "based dynamics}, 4 {XPBD distance constraints}' group=Simulation");
    TwType enumType3 = TwDefineEnum("BendingMethodType", nullptr, 0);
    TwAddVarCB(MiniGL::getTweakBar(), "BendingMethod", enumType3, setBendingMethod, getBendingMethod, &bendingMethod,
               " label='Bending method' enum='0 {None}, 1 {Dihedral angle}, 2 {Isometric bending}, 3 {XPBD isometric "
               "bending}' group=Bending");
    TwAddVarCB(MiniGL::getTweakBar(), "BendingStiffness", TW_TYPE_REAL, setBendingStiffness, getBendingStiffness, model,
               " label='Bending stiffness'  min=0.0 step=0.1 precision=4 group='Bending' ");
    TwAddVarCB(MiniGL::getTweakBar(), "DistanceStiffness", TW_TYPE_REAL, setDistanceStiffness, getDistanceStiffness,
               model, " label='Distance constraint stiffness'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "xxStiffness", TW_TYPE_REAL, setXXStiffness, getXXStiffness, model,
               " label='xx stiffness'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "yyStiffness", TW_TYPE_REAL, setYYStiffness, getYYStiffness, model,
               " label='yy stiffness'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "xyStiffness", TW_TYPE_REAL, setXYStiffness, getXYStiffness, model,
               " label='xy stiffness'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "xyPoissonRatio", TW_TYPE_REAL, setXYPoissonRatio, getXYPoissonRatio, model,
               " label='xy Poisson ratio'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "yxPoissonRatio", TW_TYPE_REAL, setYXPoissonRatio, getYXPoissonRatio, model,
               " label='yx Poisson ratio'  min=0.0 step=0.1 precision=4 group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "normalizeStretch", TW_TYPE_BOOL32, setNormalizeStretch, getNormalizeStretch,
               model, " label='Normalize stretch' group='Cloth' ");
    TwAddVarCB(MiniGL::getTweakBar(), "normalizeShear", TW_TYPE_BOOL32, setNormalizeShear, getNormalizeShear, model,
               " label='Normalize shear' group='Cloth' ");

    MiniGL::mainLoop();

    base->cleanup();

    Timing::printAverageTimes();
    Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

void reset() {
    Timing::printAverageTimes();
    Timing::reset();

    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    cd.cleanup();

    buildModel();
}

void timeStep() {
    const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
    if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime())) base->setValue(DemoBase::PAUSE, true);

    if (base->getValue<bool>(DemoBase::PAUSE)) return;

    // Simulation code
    SimulationModel* model = Simulation::getCurrent()->getModel();
    const auto numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
    for (unsigned int i = 0; i < numSteps; i++) {
        START_TIMING("SimStep")
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG
    }

    for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
        model->getTriangleModels()[i]->updateMeshNormals(model->getParticles());
}

void loadObj(const std::string& filename, VertexData& vd, utility::IndexedFaceMesh& mesh, const Vector3r& scale) {
    std::vector<utility::OBJLoader::Vec3f> x;
    std::vector<utility::OBJLoader::Vec3f> normals;
    std::vector<utility::OBJLoader::Vec2f> texCoords;
    std::vector<utility::MeshFaceIndices> faces;
    utility::OBJLoader::Vec3f s = {(float)scale[0], (float)scale[1], (float)scale[2]};
    utility::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

    mesh.release();
    const unsigned int nPoints = (unsigned int)x.size();
    const unsigned int nFaces = (unsigned int)faces.size();
    const unsigned int nTexCoords = (unsigned int)texCoords.size();
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

void buildModel() {
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

    createMesh();

    // create static rigid body
    string fileName = utility::FileSystem::normalizePath(base->getExePath() + "/resources/models/cube.obj");
    utility::IndexedFaceMesh mesh;
    VertexData vd;
    loadObj(fileName, vd, mesh, Vector3r::Ones());
    mesh.setFlatShading(true);

    string fileNameTorus = utility::FileSystem::normalizePath(base->getExePath() + "/resources/models/torus.obj");
    utility::IndexedFaceMesh meshTorus;
    VertexData vdTorus;
    loadObj(fileNameTorus, vdTorus, meshTorus, Vector3r::Ones());

    SimulationModel* model = Simulation::getCurrent()->getModel();
    SimulationModel::RigidBodyVector& rb = model->getRigidBodies();
    rb.resize(2);

    // floor
    rb[0] = new RigidBody();
    rb[0]->initBody(1.0, Vector3r(0.0, -2.5, 0.0), Quaternionr(1.0, 0.0, 0.0, 0.0), vd, mesh,
                    Vector3r(100.0, 1.0, 100.0));
    rb[0]->setMass(0.0);

    // torus
    rb[1] = new RigidBody();
    rb[1]->initBody(1.0, Vector3r(0.0, 1.5, 0.0), Quaternionr(1.0, 0.0, 0.0, 0.0), vdTorus, meshTorus,
                    Vector3r(2.0, 2.0, 2.0));
    rb[1]->setMass(0.0);
    rb[1]->setFrictionCoeff(static_cast<Real>(0.1));

    Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);
    cd.setTolerance(static_cast<Real>(0.05));

    const std::vector<Vector3r>& vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
    const unsigned int nVert1 = static_cast<unsigned int>(vertices1.size());
    cd.addCollisionBox(0, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices1.data(), nVert1,
                       Vector3r(100.0, 1.0, 100.0));

    const std::vector<Vector3r>& vertices2 = rb[1]->getGeometry().getVertexDataLocal().getVertices();
    const unsigned int nVert2 = static_cast<unsigned int>(vertices2.size());
    cd.addCollisionTorus(1, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, vertices2.data(), nVert2,
                         Vector2r(2.0, 1.0));

    SimulationModel::TriangleModelVector& tm = model->getTriangleModels();
    ParticleData& pd = model->getParticles();
    for (unsigned int i = 0; i < tm.size(); i++) {
        const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
        unsigned int offset = tm[i]->getIndexOffset();
        tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
        cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TriangleModelCollisionObjectType,
                                             &pd.getPosition(offset), nVert, true);
    }
}

void render() { base->render(); }

/** Create a particle model mesh
 */
void createMesh() {
    SimulationModel* model = Simulation::getCurrent()->getModel();
    model->addRegularTriangleModel(nCols, nRows, Vector3r(-5, 4, -5),
                                   AngleAxisr(M_PI * 0.5, Vector3r(1, 0, 0)).matrix(), Vector2r(width, height));

    // init constraints
    for (unsigned int cm = 0; cm < model->getTriangleModels().size(); cm++) {
        distanceStiffness = 1.0;
        if (simulationMethod == 4) distanceStiffness = 100000;
        model->addClothConstraints(model->getTriangleModels()[cm], simulationMethod, distanceStiffness, xxStiffness,
                                   yyStiffness, xyStiffness, xyPoissonRatio, yxPoissonRatio, normalizeStretch,
                                   normalizeShear);

        bendingStiffness = 0.01;
        if (bendingMethod == 3) bendingStiffness = 100.0;
        model->addBendingConstraints(model->getTriangleModels()[cm], bendingMethod, bendingStiffness);
    }

    LOGI("Number of triangles: {}", model->getTriangleModels()[0]->getParticleMesh().numFaces())
    LOGI("Number of vertices: {}", nRows * nCols)
}

void TW_CALL setBendingMethod(const void* value, void* clientData) {
    const short val = *(const short*)(value);
    *((short*)clientData) = val;
    reset();
}

void TW_CALL getBendingMethod(void* value, void* clientData) { *(short*)(value) = *((short*)clientData); }

void TW_CALL setSimulationMethod(const void* value, void* clientData) {
    const short val = *(const short*)(value);
    *((short*)clientData) = val;
    reset();
}

void TW_CALL getSimulationMethod(void* value, void* clientData) { *(short*)(value) = *((short*)clientData); }

void TW_CALL setBendingStiffness(const void* value, void* clientData) {
    bendingStiffness = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<DihedralConstraint, Real, &DihedralConstraint::m_stiffness>(bendingStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<IsometricBendingConstraint, Real, &IsometricBendingConstraint::m_stiffness>(
                    bendingStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<IsometricBendingConstraint_XPBD, Real, &IsometricBendingConstraint_XPBD::m_stiffness>(
                    bendingStiffness);
}

void TW_CALL getBendingStiffness(void* value, void* clientData) { *(Real*)(value) = bendingStiffness; }

void TW_CALL setDistanceStiffness(const void* value, void* clientData) {
    distanceStiffness = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<DistanceConstraint, Real, &DistanceConstraint::m_stiffness>(distanceStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<DistanceConstraint_XPBD, Real, &DistanceConstraint_XPBD::m_stiffness>(
                    distanceStiffness);
}

void TW_CALL getDistanceStiffness(void* value, void* clientData) { *(Real*)(value) = distanceStiffness; }

void TW_CALL setXXStiffness(const void* value, void* clientData) {
    xxStiffness = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<FEMTriangleConstraint, Real, &FEMTriangleConstraint::m_xxStiffness>(xxStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<StrainTriangleConstraint, Real, &StrainTriangleConstraint::m_xxStiffness>(xxStiffness);
}

void TW_CALL getXXStiffness(void* value, void* clientData) { *(Real*)(value) = xxStiffness; }

void TW_CALL getYYStiffness(void* value, void* clientData) { *(Real*)(value) = yyStiffness; }

void TW_CALL setYYStiffness(const void* value, void* clientData) {
    yyStiffness = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<FEMTriangleConstraint, Real, &FEMTriangleConstraint::m_yyStiffness>(yyStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<StrainTriangleConstraint, Real, &StrainTriangleConstraint::m_yyStiffness>(yyStiffness);
}

void TW_CALL getXYStiffness(void* value, void* clientData) { *(Real*)(value) = xyStiffness; }

void TW_CALL setXYStiffness(const void* value, void* clientData) {
    xyStiffness = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<FEMTriangleConstraint, Real, &FEMTriangleConstraint::m_xyStiffness>(xyStiffness);
    ((SimulationModel*)clientData)
            ->setConstraintValue<StrainTriangleConstraint, Real, &StrainTriangleConstraint::m_xyStiffness>(xyStiffness);
}

void TW_CALL getXYPoissonRatio(void* value, void* clientData) { *(Real*)(value) = xyPoissonRatio; }

void TW_CALL setXYPoissonRatio(const void* value, void* clientData) {
    xyPoissonRatio = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<FEMTriangleConstraint, Real, &FEMTriangleConstraint::m_xyPoissonRatio>(xyPoissonRatio);
}

void TW_CALL getYXPoissonRatio(void* value, void* clientData) { *(Real*)(value) = yxPoissonRatio; }

void TW_CALL setYXPoissonRatio(const void* value, void* clientData) {
    yxPoissonRatio = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<FEMTriangleConstraint, Real, &FEMTriangleConstraint::m_yxPoissonRatio>(yxPoissonRatio);
}

void TW_CALL getNormalizeStretch(void* value, void* clientData) { *(bool*)(value) = normalizeStretch; }

void TW_CALL setNormalizeStretch(const void* value, void* clientData) {
    normalizeStretch = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<StrainTriangleConstraint, bool, &StrainTriangleConstraint::m_normalizeStretch>(
                    normalizeStretch);
}

void TW_CALL getNormalizeShear(void* value, void* clientData) { *(bool*)(value) = normalizeShear; }

void TW_CALL setNormalizeShear(const void* value, void* clientData) {
    normalizeShear = *(const Real*)(value);
    ((SimulationModel*)clientData)
            ->setConstraintValue<StrainTriangleConstraint, bool, &StrainTriangleConstraint::m_normalizeShear>(
                    normalizeShear);
}