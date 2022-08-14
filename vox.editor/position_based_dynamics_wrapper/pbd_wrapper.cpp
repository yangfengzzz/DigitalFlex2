//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/position_based_dynamics_wrapper/pbd_wrapper.h"

#include <cmath>
#include <iostream>

#include "vox.base/geometry/triangle_mesh_distance.h"
#include "vox.base/mesh/triangle_mesh.h"
#include "vox.pbd/common.h"
#include "vox.pbd/file_system.h"
#include "vox.pbd/obj_loader.h"
#include "vox.pbd/scene_loader.h"
#include "vox.pbd/simulation.h"
#include "vox.pbd/tet_gen_loader.h"
#include "vox.pbd/time_manager.h"
#include "vox.pbd/time_step.h"
#include "vox.pbd/time_step_controller.h"
#include "vox.pbd/timing.h"

using namespace Eigen;
using namespace std;

namespace vox {

PBDWrapper::PBDWrapper() {
    m_clothSimulationMethod = 2;
    m_solidSimulationMethod = 2;
    m_bendingMethod = 2;
    m_sceneName = "";
    m_sceneFileName = "";
    m_dampingCoeff = 0.0;
    m_timeStep = new TimeStepController();
    m_timeStep->init();
    m_model.init();
    Simulation::getCurrent()->setModel(&m_model);
}

PBDWrapper::~PBDWrapper() {
    delete TimeManager::getCurrent();
    delete m_timeStep;
}

// void PBDWrapper::initGUI_TweakBar(TwBar *tweakBar)
//{
//	TwType enumType2 = TwDefineEnum("ClothSimulationMethodType", NULL, 0);
//	TwAddVarCB(tweakBar, "ClothSimulationMethod", enumType2, setClothSimulationMethod, getClothSimulationMethod,
// this, " label='Cloth sim. method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain based
// dynamics}' group=PBD"); 	TwType enumType3 = TwDefineEnum("SolidSimulationMethodType", NULL, 0);
// TwAddVarCB(tweakBar, "SolidSimulationMethod", enumType3, setSolidSimulationMethod, getSolidSimulationMethod, this,
//" label='Solid sim. method' enum='0 {None}, 1 {Volume constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics (no
// inversion handling)}, 4 {Shape matching (no inversion handling)}' group=PBD"); 	TwAddVarCB(tweakBar,
// "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, &m_model, " label='Stiffness'  min=0.0 step=0.1 precision=4
// group='Distance constraints' "); 	TwAddVarCB(tweakBar, "XXStiffness", TW_TYPE_REAL, setXXStiffness,
// getXXStiffness,
//&m_model, " label='Stiffness
// XX'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' "); 	TwAddVarCB(tweakBar, "YYStiffness",
// TW_TYPE_REAL, setYYStiffness, getYYStiffness, &m_model, " label='Stiffness YY'  min=0.0 step=0.1 precision=4
// group='Strain based
// dynamics' "); 	TwAddVarCB(tweakBar, "XYStiffness", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &m_model, "
// label='Stiffness XY'  min=0.0 step=0.1 precision=4 group='Strain based dynamics' "); 	TwAddVarCB(tweakBar,
//"XXStiffnessFEM", TW_TYPE_REAL, setXXStiffness, getXXStiffness, &m_model, " label='Youngs modulus XX'  min=0.0
// step=0.1 precision=4 group='FEM based PBD' "); 	TwAddVarCB(tweakBar, "YYStiffnessFEM", TW_TYPE_REAL,
// setYYStiffness, getYYStiffness, &m_model, " label='Youngs modulus YY'  min=0.0 step=0.1 precision=4 group='FEM based
// PBD' "); 	TwAddVarCB(tweakBar, "XYStiffnessFEM", TW_TYPE_REAL, setXYStiffness, getXYStiffness, &m_model, "
// label='Youngs modulus XY'  min=0.0 step=0.1 precision=4 group='FEM based PBD' "); 	TwAddVarCB(tweakBar,
// "XYPoissonRatioFEM", TW_TYPE_REAL, setXYPoissonRatio, getXYPoissonRatio, &m_model, " label='Poisson ratio XY' min=0.0
// step=0.1 precision=4 group='FEM based PBD' "); 	TwAddVarCB(tweakBar, "YXPoissonRatioFEM", TW_TYPE_REAL,
// setYXPoissonRatio, getYXPoissonRatio, &m_model, " label='Poisson ratio YX'  min=0.0 step=0.1 precision=4 group='FEM
// based PBD' "); 	TwAddVarCB(tweakBar, "NormalizeStretch", TW_TYPE_BOOL32, setNormalizeStretch,
// getNormalizeStretch,
//&m_model, " label='Normalize stretch' group='Strain based dynamics' "); 	TwAddVarCB(tweakBar, "NormalizeShear",
// TW_TYPE_BOOL32, setNormalizeShear, getNormalizeShear, &m_model, " label='Normalize shear' group='Strain based
// dynamics' "); 	TwType enumType4 = TwDefineEnum("BendingMethodType", NULL, 0); 	TwAddVarCB(tweakBar,
// "BendingMethod", enumType4, setBendingMethod, getBendingMethod, this, " label='Bending method' enum='0 {None}, 1
// {Dihedral angle}, 2 {Isometric bending}' group=Bending"); 	TwAddVarCB(tweakBar, "BendingStiffness", TW_TYPE_REAL,
// setBendingStiffness, getBendingStiffness, &m_model, " label='Bending stiffness'  min=0.0 step=0.01 precision=4
// group=Bending ");
// }
void PBDWrapper::reset() {
    m_model.reset();
    m_timeStep->reset();
}

void PBDWrapper::timeStep() {
    ParticleData &pd = m_model.getParticles();
    SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
    TimeManager::getCurrent()->setTimeStepSize(TimeManager::getCurrent()->getTimeStepSize());
    TimeManager::getCurrent()->setTime(TimeManager::getCurrent()->getTime());

    m_timeStep->step(m_model);

    for (unsigned int i = 0; i < pd.size(); i++) {
        pd.getVelocity(i) *= (static_cast<Real>(1.0) - m_dampingCoeff);
    }
    for (unsigned int i = 0; i < rb.size(); i++) {
        rb[i]->getVelocity() *= (static_cast<Real>(1.0) - m_dampingCoeff);
        rb[i]->getAngularVelocity() *= (static_cast<Real>(1.0) - m_dampingCoeff);
    }
}

void PBDWrapper::updateVisModels() {
    ParticleData &pd = m_model.getParticles();

    // Update visualization models
    for (unsigned int i = 0; i < m_model.getTetModels().size(); i++) {
        m_model.getTetModels()[i]->updateMeshNormals(pd);
        m_model.getTetModels()[i]->updateVisMesh(pd);
    }
    for (unsigned int i = 0; i < m_model.getTriangleModels().size(); i++) {
        m_model.getTriangleModels()[i]->updateMeshNormals(pd);
    }
}

void PBDWrapper::loadObj(const std::string &filename,
                         VertexData &vd,
                         utility::IndexedFaceMesh &mesh,
                         const Vector3r &scale) {
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

void PBDWrapper::readScene(const std::string &sceneFileName, const std::vector<RBData> &additionalRigidBodies) {
    m_sceneFileName = sceneFileName;

    SimulationModel::RigidBodyVector &rb = m_model.getRigidBodies();
    SimulationModel::TriangleModelVector &triModels = m_model.getTriangleModels();
    SimulationModel::TetModelVector &tetModels = m_model.getTetModels();
    SimulationModel::ConstraintVector &constraints = m_model.getConstraints();

    utility::SceneLoader::SceneData data;
    utility::SceneLoader loader;
    loader.readScene(sceneFileName, data);
    LOGI("Scene: {}", sceneFileName)

    std::string basePath = utility::FileSystem::getFilePath(sceneFileName);

    m_sceneName = data.m_sceneName;

    Simulation *sim = Simulation::getCurrent();
    sim->setVecValue<Real>(Simulation::GRAVITATION, &data.m_gravity[0]);
    m_timeStep->setValue(TimeStepController::MAX_ITERATIONS, data.m_maxIter);
    m_timeStep->setValue(TimeStepController::MAX_ITERATIONS_V, data.m_maxIterVel);
    m_timeStep->setValue(TimeStepController::VELOCITY_UPDATE_METHOD, data.m_velocityUpdateMethod);
    if (data.m_triangleModelSimulationMethod != -1) m_clothSimulationMethod = data.m_triangleModelSimulationMethod;
    if (data.m_tetModelSimulationMethod != -1) m_solidSimulationMethod = data.m_tetModelSimulationMethod;
    if (data.m_triangleModelBendingMethod != -1) m_bendingMethod = data.m_triangleModelBendingMethod;
    m_cd.setTolerance(data.m_contactTolerance);
    m_model.setContactStiffnessRigidBody(data.m_contactStiffnessRigidBody);
    m_model.setContactStiffnessParticleRigidBody(data.m_contactStiffnessParticleRigidBody);

    m_model.setValue(SimulationModel::CLOTH_STIFFNESS, data.m_cloth_stiffness);
    m_model.setValue(SimulationModel::CLOTH_BENDING_STIFFNESS, data.m_cloth_bendingStiffness);
    m_model.setValue(SimulationModel::CLOTH_STIFFNESS_XX, data.m_cloth_xxStiffness);
    m_model.setValue(SimulationModel::CLOTH_STIFFNESS_YY, data.m_cloth_yyStiffness);
    m_model.setValue(SimulationModel::CLOTH_STIFFNESS_XY, data.m_cloth_xyStiffness);
    m_model.setValue(SimulationModel::CLOTH_POISSON_RATIO_XY, data.m_cloth_xyPoissonRatio);
    m_model.setValue(SimulationModel::CLOTH_POISSON_RATIO_YX, data.m_cloth_yxPoissonRatio);
    m_model.setValue(SimulationModel::CLOTH_NORMALIZE_STRETCH, data.m_cloth_normalizeStretch);
    m_model.setValue(SimulationModel::CLOTH_NORMALIZE_SHEAR, data.m_cloth_normalizeShear);

    //////////////////////////////////////////////////////////////////////////
    // rigid bodies
    //////////////////////////////////////////////////////////////////////////

    // add additional bodies
    for (auto i = 0; i < additionalRigidBodies.size(); i++) {
        utility::SceneLoader::RigidBodyData rbd;
        rbd.m_modelFile = additionalRigidBodies[i].objFile;
        rbd.m_collisionObjectFileName = "";
        rbd.m_x = additionalRigidBodies[i].x;
        rbd.m_q = Quaternionr(additionalRigidBodies[i].R);
        rbd.m_scale = additionalRigidBodies[i].scale;
        rbd.m_collisionObjectScale = additionalRigidBodies[i].scale;
        rbd.m_collisionObjectType = additionalRigidBodies[i].collisionType;
        rbd.m_frictionCoeff = additionalRigidBodies[i].friction;
        rbd.m_restitutionCoeff = additionalRigidBodies[i].restitution;
        rbd.m_density = 1000.0;
        rbd.m_invertSDF = false;
        rbd.m_isDynamic = false;
        rbd.m_omega.setZero();
        rbd.m_v.setZero();
        rbd.m_resolutionSDF = Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>(10, 10, 10);
        rbd.m_thicknessSDF = 0.0;
        rbd.m_testMesh = false;
        data.m_rigidBodyData.push_back(rbd);
    }

    // map file names to loaded geometry to prevent multiple imports of same files
    std::map<std::string, pair<VertexData, utility::IndexedFaceMesh>> objFiles;
    std::map<std::string, CubicSDFCollisionDetection::GridPtr> distanceFields;
    for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++) {
        utility::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

        // Check if already loaded
        rbd.m_modelFile = utility::FileSystem::normalizePath(rbd.m_modelFile);
        if (objFiles.find(rbd.m_modelFile) == objFiles.end()) {
            utility::IndexedFaceMesh mesh;
            VertexData vd;
            loadObj(rbd.m_modelFile, vd, mesh, Vector3r::Ones());
            objFiles[rbd.m_modelFile] = {vd, mesh};
        }

        const std::string basePath = utility::FileSystem::getFilePath(sceneFileName);
        const string cachePath = basePath + "/Cache";
        const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) + "_" +
                              to_string(rbd.m_resolutionSDF[2]);
        const std::string modelFileName = utility::FileSystem::getFileNameWithExt(rbd.m_modelFile);
        const string sdfFileName =
                utility::FileSystem::normalizePath(cachePath + "/" + modelFileName + "_" + resStr + ".csdf");

        std::string sdfKey = rbd.m_collisionObjectFileName;
        if (sdfKey == "") {
            sdfKey = sdfFileName;
        }
        if (distanceFields.find(sdfKey) == distanceFields.end()) {
            // Generate SDF
            if (rbd.m_collisionObjectType == utility::SceneLoader::SDF) {
                if (rbd.m_collisionObjectFileName == "") {
                    std::string md5FileName =
                            utility::FileSystem::normalizePath(cachePath + "/" + modelFileName + ".md5");
                    string md5Str = utility::FileSystem::getFileMD5(rbd.m_modelFile);
                    bool md5 = false;
                    if (utility::FileSystem::fileExists(md5FileName))
                        md5 = utility::FileSystem::checkMD5(md5Str, md5FileName);

                    // check MD5 if cache file is available
                    const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) +
                                          "_" + to_string(rbd.m_resolutionSDF[2]);
                    const string sdfFileName = utility::FileSystem::normalizePath(cachePath + "/" + modelFileName +
                                                                                  "_" + resStr + ".csdf");
                    bool foundCacheFile = utility::FileSystem::fileExists(sdfFileName);

                    if (foundCacheFile && md5) {
                        LOGI("Load cached SDF: {}", sdfFileName)
                        distanceFields[sdfFileName] = std::make_shared<CubicSDFCollisionDetection::Grid>(sdfFileName);
                    } else {
                        VertexData &vd = objFiles[rbd.m_modelFile].first;
                        utility::IndexedFaceMesh &mesh = objFiles[rbd.m_modelFile].second;

                        std::vector<unsigned int> &faces = mesh.getFaces();
                        const unsigned int nFaces = mesh.numFaces();

#ifdef USE_DOUBLE
                        TriangleMesh sdfMesh(&vd.getPosition(0)[0], faces.data(), vd.size(), nFaces);
#else
                        // if type is float, copy vector to double vector
                        std::vector<double> doubleVec;
                        doubleVec.resize(3 * vd.size());
                        for (unsigned int i = 0; i < vd.size(); i++)
                            for (unsigned int j = 0; j < 3; j++) doubleVec[3 * i + j] = vd.getPosition(i)[j];
                        TriangleMesh sdfMesh(&doubleVec[0], faces.data(), vd.size(), nFaces);
#endif
                        TriangleMeshDistance md(sdfMesh);
                        Eigen::AlignedBox3d domain;
                        for (auto const &x : sdfMesh.vertices()) {
                            domain.extend(x);
                        }
                        domain.max() += 0.1 * Eigen::Vector3d::Ones();
                        domain.min() -= 0.1 * Eigen::Vector3d::Ones();

                        LOGI("Set SDF resolution: {}, {}, {}", rbd.m_resolutionSDF[0], rbd.m_resolutionSDF[1],
                             rbd.m_resolutionSDF[2])
                        distanceFields[sdfFileName] = std::make_shared<CubicSDFCollisionDetection::Grid>(
                                domain, std::array<unsigned int, 3>({rbd.m_resolutionSDF[0], rbd.m_resolutionSDF[1],
                                                                     rbd.m_resolutionSDF[2]}));
                        auto func = DiscreteGrid::ContinuousFunction{};
                        func = [&md](Eigen::Vector3d const &xi) { return md.signed_distance(xi); };
                        LOGI("Generate SDF for {}", rbd.m_modelFile)
                        distanceFields[sdfFileName]->addFunction(func, true);
                        if (utility::FileSystem::makeDir(cachePath) == 0) {
                            LOGI("Save SDF: {}", sdfFileName)
                            distanceFields[sdfFileName]->save(sdfFileName);
                            utility::FileSystem::writeMD5File(rbd.m_modelFile, md5FileName);
                        }
                    }
                } else {
                    std::string fileName = rbd.m_collisionObjectFileName;
                    if (utility::FileSystem::isRelativePath(fileName)) {
                        fileName = utility::FileSystem::normalizePath(basePath + "/" + fileName);
                    }
                    LOGI("Load SDF: {}", fileName)
                    distanceFields[rbd.m_collisionObjectFileName] =
                            std::make_shared<CubicSDFCollisionDetection::Grid>(fileName);
                }
            }
        }
    }

    for (unsigned int i = 0; i < data.m_tetModelData.size(); i++) {
        const utility::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

        // Check if already loaded
        if ((tmd.m_modelFileVis != "") && (objFiles.find(tmd.m_modelFileVis) == objFiles.end())) {
            utility::IndexedFaceMesh mesh;
            VertexData vd;
            loadObj(utility::FileSystem::normalizePath(tmd.m_modelFileVis), vd, mesh, Vector3r::Ones());
            objFiles[tmd.m_modelFileVis] = {vd, mesh};
        }
    }

    rb.resize(data.m_rigidBodyData.size());
    std::map<unsigned int, unsigned int> id_index;
    for (unsigned int i = 0; i < data.m_rigidBodyData.size(); i++) {
        const utility::SceneLoader::RigidBodyData &rbd = data.m_rigidBodyData[i];

        if (objFiles.find(rbd.m_modelFile) == objFiles.end()) continue;

        id_index[rbd.m_id] = i;

        VertexData &vd = objFiles[rbd.m_modelFile].first;
        utility::IndexedFaceMesh &mesh = objFiles[rbd.m_modelFile].second;

        rb[i] = new RigidBody();

        rb[i]->initBody(rbd.m_density, rbd.m_x, rbd.m_q, vd, mesh, rbd.m_scale);

        if (!rbd.m_isDynamic) {
            rb[i]->setMass(0.0);
            rb[i]->setVelocity(Vector3r::Zero());
            rb[i]->setVelocity0(Vector3r::Zero());
            rb[i]->setAngularVelocity(Vector3r::Zero());
            rb[i]->setAngularVelocity0(Vector3r::Zero());
        } else {
            rb[i]->setVelocity(rbd.m_v);
            rb[i]->setVelocity0(rbd.m_v);
            rb[i]->setAngularVelocity(rbd.m_omega);
            rb[i]->setAngularVelocity0(rbd.m_omega);
        }
        rb[i]->setRestitutionCoeff(rbd.m_restitutionCoeff);
        rb[i]->setFrictionCoeff(rbd.m_frictionCoeff);

        const std::vector<Vector3r> &vertices = rb[i]->getGeometry().getVertexDataLocal().getVertices();
        const unsigned int nVert = static_cast<unsigned int>(vertices.size());

        switch (rbd.m_collisionObjectType) {
            case utility::SceneLoader::No_Collision_Object:
                break;
            case utility::SceneLoader::Sphere:
                m_cd.addCollisionSphere(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                        &(vertices)[0], nVert, rbd.m_collisionObjectScale[0], rbd.m_testMesh,
                                        rbd.m_invertSDF);
                break;
            case utility::SceneLoader::Box:
                m_cd.addCollisionBox(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                     &(vertices)[0], nVert, rbd.m_collisionObjectScale, rbd.m_testMesh,
                                     rbd.m_invertSDF);
                break;
            case utility::SceneLoader::Cylinder:
                m_cd.addCollisionCylinder(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                          &(vertices)[0], nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh,
                                          rbd.m_invertSDF);
                break;
            case utility::SceneLoader::Torus:
                m_cd.addCollisionTorus(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                       &(vertices)[0], nVert, rbd.m_collisionObjectScale.head<2>(), rbd.m_testMesh,
                                       rbd.m_invertSDF);
                break;
            case utility::SceneLoader::HollowSphere:
                m_cd.addCollisionHollowSphere(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                              &(vertices)[0], nVert, rbd.m_collisionObjectScale[0], rbd.m_thicknessSDF,
                                              rbd.m_testMesh, rbd.m_invertSDF);
                break;
            case utility::SceneLoader::HollowBox:
                m_cd.addCollisionHollowBox(i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                           &(vertices)[0], nVert, rbd.m_collisionObjectScale, rbd.m_thicknessSDF,
                                           rbd.m_testMesh, rbd.m_invertSDF);
                break;
            case utility::SceneLoader::SDF: {
                if (rbd.m_collisionObjectFileName == "") {
                    const std::string basePath = utility::FileSystem::getFilePath(sceneFileName);
                    const string cachePath = basePath + "/Cache";
                    const string resStr = to_string(rbd.m_resolutionSDF[0]) + "_" + to_string(rbd.m_resolutionSDF[1]) +
                                          "_" + to_string(rbd.m_resolutionSDF[2]);
                    const std::string modelFileName = utility::FileSystem::getFileNameWithExt(rbd.m_modelFile);
                    const string sdfFileName = utility::FileSystem::normalizePath(cachePath + "/" + modelFileName +
                                                                                  "_" + resStr + ".csdf");
                    m_cd.addCubicSDFCollisionObject(
                            i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(vertices)[0], nVert,
                            distanceFields[sdfFileName], rbd.m_collisionObjectScale, rbd.m_testMesh, rbd.m_invertSDF);
                } else
                    m_cd.addCubicSDFCollisionObject(
                            i, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(vertices)[0], nVert,
                            distanceFields[rbd.m_collisionObjectFileName], rbd.m_collisionObjectScale, rbd.m_testMesh,
                            rbd.m_invertSDF);
                break;
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // triangle models
    //////////////////////////////////////////////////////////////////////////

    // map file names to loaded geometry to prevent multiple imports of same files
    std::map<std::string, pair<VertexData, utility::IndexedFaceMesh>> triFiles;
    for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++) {
        const utility::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

        // Check if already loaded
        if (triFiles.find(tmd.m_modelFile) == triFiles.end()) {
            utility::IndexedFaceMesh mesh;
            VertexData vd;
            loadObj(utility::FileSystem::normalizePath(tmd.m_modelFile), vd, mesh, Vector3r::Ones());
            triFiles[tmd.m_modelFile] = {vd, mesh};
        }
    }

    triModels.reserve(data.m_triangleModelData.size());
    std::map<unsigned int, unsigned int> tm_id_index;
    for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++) {
        const utility::SceneLoader::TriangleModelData &tmd = data.m_triangleModelData[i];

        if (triFiles.find(tmd.m_modelFile) == triFiles.end()) continue;

        tm_id_index[tmd.m_id] = i;

        VertexData vd = triFiles[tmd.m_modelFile].first;
        utility::IndexedFaceMesh &mesh = triFiles[tmd.m_modelFile].second;

        const Matrix3r R = tmd.m_q.matrix();
        for (unsigned int j = 0; j < vd.size(); j++) {
            vd.getPosition(j) = R * (vd.getPosition(j).cwiseProduct(tmd.m_scale)) + tmd.m_x;
        }

        m_model.addTriangleModel(vd.size(), mesh.numFaces(), &vd.getPosition(0), mesh.getFaces().data(),
                                 mesh.getUVIndices(), mesh.getUVs());

        TriangleModel *tm = triModels[triModels.size() - 1];
        ParticleData &pd = m_model.getParticles();
        unsigned int offset = tm->getIndexOffset();

        for (unsigned int j = 0; j < tmd.m_staticParticles.size(); j++) {
            const unsigned int index = tmd.m_staticParticles[j] + offset;
            pd.setMass(index, 0.0);
        }

        tm->setRestitutionCoeff(tmd.m_restitutionCoeff);
        tm->setFrictionCoeff(tmd.m_frictionCoeff);
    }

    initTriangleModelConstraints();

    //////////////////////////////////////////////////////////////////////////
    // tet models
    //////////////////////////////////////////////////////////////////////////

    // map file names to loaded geometry to prevent multiple imports of same files
    std::map<pair<string, string>, pair<vector<Vector3r>, vector<unsigned int>>> tetFiles;
    for (unsigned int i = 0; i < data.m_tetModelData.size(); i++) {
        const utility::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

        // Check if already loaded
        pair<string, string> fileNames = {tmd.m_modelFileNodes, tmd.m_modelFileElements};
        if (tetFiles.find(fileNames) == tetFiles.end()) {
            vector<Vector3r> vertices;
            vector<unsigned int> tets;
            utility::TetGenLoader::loadTetgenModel(utility::FileSystem::normalizePath(tmd.m_modelFileNodes),
                                                   utility::FileSystem::normalizePath(tmd.m_modelFileElements),
                                                   vertices, tets);
            tetFiles[fileNames] = {vertices, tets};
        }
    }

    tetModels.reserve(data.m_tetModelData.size());
    std::map<unsigned int, unsigned int> tm_id_index2;
    for (unsigned int i = 0; i < data.m_tetModelData.size(); i++) {
        const utility::SceneLoader::TetModelData &tmd = data.m_tetModelData[i];

        pair<string, string> fileNames = {tmd.m_modelFileNodes, tmd.m_modelFileElements};
        auto geo = tetFiles.find(fileNames);
        if (geo == tetFiles.end()) continue;

        tm_id_index2[tmd.m_id] = i;

        vector<Vector3r> vertices = geo->second.first;
        vector<unsigned int> &tets = geo->second.second;

        const Matrix3r R = tmd.m_q.matrix();
        for (unsigned int j = 0; j < vertices.size(); j++) {
            vertices[j] = R * (vertices[j].cwiseProduct(tmd.m_scale)) + tmd.m_x;
        }

        m_model.addTetModel((unsigned int)vertices.size(), (unsigned int)tets.size() / 4, vertices.data(), tets.data());

        TetModel *tm = tetModels[tetModels.size() - 1];
        ParticleData &pd = m_model.getParticles();
        unsigned int offset = tm->getIndexOffset();

        for (unsigned int j = 0; j < tmd.m_staticParticles.size(); j++) {
            const unsigned int index = tmd.m_staticParticles[j] + offset;
            pd.setMass(index, 0.0);
        }

        // read visualization mesh
        if (tmd.m_modelFileVis != "") {
            if (objFiles.find(tmd.m_modelFileVis) != objFiles.end()) {
                utility::IndexedFaceMesh &visMesh = tm->getVisMesh();
                VertexData &vdVis = tm->getVisVertices();
                vdVis = objFiles[tmd.m_modelFileVis].first;
                visMesh = objFiles[tmd.m_modelFileVis].second;

                for (unsigned int j = 0; j < vdVis.size(); j++)
                    vdVis.getPosition(j) = R * (vdVis.getPosition(j).cwiseProduct(tmd.m_scale)) + tmd.m_x;

                tm->updateMeshNormals(pd);
                tm->attachVisMesh(pd);
                tm->updateVisMesh(pd);
            }
        }

        tm->setRestitutionCoeff(tmd.m_restitutionCoeff);
        tm->setFrictionCoeff(tmd.m_frictionCoeff);

        tm->updateMeshNormals(pd);
    }

    initTetModelConstraints();

    // init collision objects for deformable models
    ParticleData &pd = m_model.getParticles();
    for (unsigned int i = 0; i < data.m_triangleModelData.size(); i++) {
        TriangleModel *tm = triModels[i];
        unsigned int offset = tm->getIndexOffset();
        const unsigned int nVert = tm->getParticleMesh().numVertices();
        m_cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TriangleModelCollisionObjectType,
                                               &pd.getPosition(offset), nVert, true);
    }
    for (unsigned int i = 0; i < data.m_tetModelData.size(); i++) {
        TetModel *tm = tetModels[i];
        unsigned int offset = tm->getIndexOffset();
        const unsigned int nVert = tm->getParticleMesh().numVertices();
        m_cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TetModelCollisionObjectType,
                                               &pd.getPosition(offset), nVert, true);
    }

    //////////////////////////////////////////////////////////////////////////
    // joints
    //////////////////////////////////////////////////////////////////////////

    for (unsigned int i = 0; i < data.m_ballJointData.size(); i++) {
        const utility::SceneLoader::BallJointData &jd = data.m_ballJointData[i];
        m_model.addBallJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position);
    }

    for (unsigned int i = 0; i < data.m_ballOnLineJointData.size(); i++) {
        const utility::SceneLoader::BallOnLineJointData &jd = data.m_ballOnLineJointData[i];
        m_model.addBallOnLineJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
    }

    for (unsigned int i = 0; i < data.m_hingeJointData.size(); i++) {
        const utility::SceneLoader::HingeJointData &jd = data.m_hingeJointData[i];
        m_model.addHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis);
    }

    for (unsigned int i = 0; i < data.m_universalJointData.size(); i++) {
        const utility::SceneLoader::UniversalJointData &jd = data.m_universalJointData[i];
        m_model.addUniversalJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position, jd.m_axis[0],
                                  jd.m_axis[1]);
    }

    for (unsigned int i = 0; i < data.m_sliderJointData.size(); i++) {
        const utility::SceneLoader::SliderJointData &jd = data.m_sliderJointData[i];
        m_model.addSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
    }

    for (unsigned int i = 0; i < data.m_rigidBodyParticleBallJointData.size(); i++) {
        const utility::SceneLoader::RigidBodyParticleBallJointData &jd = data.m_rigidBodyParticleBallJointData[i];
        m_model.addRigidBodyParticleBallJoint(id_index[jd.m_bodyID[0]], jd.m_bodyID[1]);
    }

    for (unsigned int i = 0; i < data.m_targetAngleMotorHingeJointData.size(); i++) {
        const utility::SceneLoader::TargetAngleMotorHingeJointData &jd = data.m_targetAngleMotorHingeJointData[i];
        m_model.addTargetAngleMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position,
                                              jd.m_axis);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTarget(jd.m_target);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
        ((MotorJoint *)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
    }

    for (unsigned int i = 0; i < data.m_targetVelocityMotorHingeJointData.size(); i++) {
        const utility::SceneLoader::TargetVelocityMotorHingeJointData &jd = data.m_targetVelocityMotorHingeJointData[i];
        m_model.addTargetVelocityMotorHingeJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position,
                                                 jd.m_axis);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTarget(jd.m_target);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
        ((MotorJoint *)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
    }

    for (unsigned int i = 0; i < data.m_targetPositionMotorSliderJointData.size(); i++) {
        const utility::SceneLoader::TargetPositionMotorSliderJointData &jd =
                data.m_targetPositionMotorSliderJointData[i];
        m_model.addTargetPositionMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTarget(jd.m_target);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
        ((MotorJoint *)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
    }

    for (unsigned int i = 0; i < data.m_targetVelocityMotorSliderJointData.size(); i++) {
        const utility::SceneLoader::TargetVelocityMotorSliderJointData &jd =
                data.m_targetVelocityMotorSliderJointData[i];
        m_model.addTargetVelocityMotorSliderJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTarget(jd.m_target);
        ((MotorJoint *)constraints[constraints.size() - 1])->setTargetSequence(jd.m_targetSequence);
        ((MotorJoint *)constraints[constraints.size() - 1])->setRepeatSequence(jd.m_repeat);
    }

    for (unsigned int i = 0; i < data.m_damperJointData.size(); i++) {
        const utility::SceneLoader::DamperJointData &jd = data.m_damperJointData[i];
        m_model.addDamperJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_axis, jd.m_stiffness);
    }

    for (unsigned int i = 0; i < data.m_rigidBodySpringData.size(); i++) {
        const utility::SceneLoader::RigidBodySpringData &jd = data.m_rigidBodySpringData[i];
        m_model.addRigidBodySpring(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position1, jd.m_position2,
                                   jd.m_stiffness);
    }

    for (unsigned int i = 0; i < data.m_distanceJointData.size(); i++) {
        const utility::SceneLoader::DistanceJointData &jd = data.m_distanceJointData[i];
        m_model.addDistanceJoint(id_index[jd.m_bodyID[0]], id_index[jd.m_bodyID[1]], jd.m_position1, jd.m_position2);
    }

    m_cd.updateAABBs(m_model);
}

void PBDWrapper::initModel(const Real timeStepSize) {
    TimeManager::getCurrent()->setTimeStepSize(timeStepSize);

    m_timeStep->setCollisionDetection(m_model, &m_cd);
}

TimeStepController &PBDWrapper::getTimeStepController() { return *m_timeStep; }

void PBDWrapper::initTriangleModelConstraints() {
    // init constraints
    for (unsigned int cm = 0; cm < m_model.getTriangleModels().size(); cm++) {
        const unsigned int offset = m_model.getTriangleModels()[cm]->getIndexOffset();
        if (m_clothSimulationMethod == 1) {
            const unsigned int nEdges = m_model.getTriangleModels()[cm]->getParticleMesh().numEdges();
            const utility::IndexedFaceMesh::Edge *edges =
                    m_model.getTriangleModels()[cm]->getParticleMesh().getEdges().data();
            for (unsigned int i = 0; i < nEdges; i++) {
                const unsigned int v1 = edges[i].m_vert[0] + offset;
                const unsigned int v2 = edges[i].m_vert[1] + offset;

                m_model.addDistanceConstraint(v1, v2);
            }
        } else if (m_clothSimulationMethod == 2) {
            TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
            const unsigned int *tris = mesh.getFaces().data();
            const unsigned int nFaces = mesh.numFaces();
            for (unsigned int i = 0; i < nFaces; i++) {
                const unsigned int v1 = tris[3 * i] + offset;
                const unsigned int v2 = tris[3 * i + 1] + offset;
                const unsigned int v3 = tris[3 * i + 2] + offset;
                m_model.addFEMTriangleConstraint(v1, v2, v3);
            }
        } else if (m_clothSimulationMethod == 3) {
            TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
            const unsigned int *tris = mesh.getFaces().data();
            const unsigned int nFaces = mesh.numFaces();
            for (unsigned int i = 0; i < nFaces; i++) {
                const unsigned int v1 = tris[3 * i] + offset;
                const unsigned int v2 = tris[3 * i + 1] + offset;
                const unsigned int v3 = tris[3 * i + 2] + offset;
                m_model.addStrainTriangleConstraint(v1, v2, v3);
            }
        }
        if (m_bendingMethod != 0) {
            TriangleModel::ParticleMesh &mesh = m_model.getTriangleModels()[cm]->getParticleMesh();
            unsigned int nEdges = mesh.numEdges();
            const TriangleModel::ParticleMesh::Edge *edges = mesh.getEdges().data();
            const unsigned int *tris = mesh.getFaces().data();
            for (unsigned int i = 0; i < nEdges; i++) {
                const int tri1 = edges[i].m_face[0];
                const int tri2 = edges[i].m_face[1];
                if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff)) {
                    // Find the triangle points which do not lie on the axis
                    const int axisPoint1 = edges[i].m_vert[0];
                    const int axisPoint2 = edges[i].m_vert[1];
                    int point1 = -1;
                    int point2 = -1;
                    for (int j = 0; j < 3; j++) {
                        if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2)) {
                            point1 = tris[3 * tri1 + j];
                            break;
                        }
                    }
                    for (int j = 0; j < 3; j++) {
                        if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2)) {
                            point2 = tris[3 * tri2 + j];
                            break;
                        }
                    }
                    if ((point1 != -1) && (point2 != -1)) {
                        const unsigned int vertex1 = point1 + offset;
                        const unsigned int vertex2 = point2 + offset;
                        const unsigned int vertex3 = edges[i].m_vert[0] + offset;
                        const unsigned int vertex4 = edges[i].m_vert[1] + offset;
                        if (m_bendingMethod == 1)
                            m_model.addDihedralConstraint(vertex1, vertex2, vertex3, vertex4);
                        else if (m_bendingMethod == 2)
                            m_model.addIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
                    }
                }
            }
        }
    }
}

void PBDWrapper::initTetModelConstraints() {
    // init constraints
    for (unsigned int cm = 0; cm < m_model.getTetModels().size(); cm++) {
        const unsigned int offset = m_model.getTetModels()[cm]->getIndexOffset();
        const unsigned int nTets = m_model.getTetModels()[cm]->getParticleMesh().numTets();
        const unsigned int *tets = m_model.getTetModels()[cm]->getParticleMesh().getTets().data();
        const utility::IndexedTetMesh::VerticesTets *vTets =
                m_model.getTetModels()[cm]->getParticleMesh().getVertexTets().data();
        if (m_solidSimulationMethod == 1) {
            const unsigned int offset = m_model.getTetModels()[cm]->getIndexOffset();
            const unsigned int nEdges = m_model.getTetModels()[cm]->getParticleMesh().numEdges();
            const utility::IndexedTetMesh::Edge *edges =
                    m_model.getTetModels()[cm]->getParticleMesh().getEdges().data();
            for (unsigned int i = 0; i < nEdges; i++) {
                const unsigned int v1 = edges[i].m_vert[0] + offset;
                const unsigned int v2 = edges[i].m_vert[1] + offset;

                m_model.addDistanceConstraint(v1, v2);
            }

            for (unsigned int i = 0; i < nTets; i++) {
                const unsigned int v1 = tets[4 * i] + offset;
                const unsigned int v2 = tets[4 * i + 1] + offset;
                const unsigned int v3 = tets[4 * i + 2] + offset;
                const unsigned int v4 = tets[4 * i + 3] + offset;

                m_model.addVolumeConstraint(v1, v2, v3, v4);
            }
        } else if (m_solidSimulationMethod == 2) {
            TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++) {
                const unsigned int v1 = tets[4 * i] + offset;
                const unsigned int v2 = tets[4 * i + 1] + offset;
                const unsigned int v3 = tets[4 * i + 2] + offset;
                const unsigned int v4 = tets[4 * i + 3] + offset;

                m_model.addFEMTetConstraint(v1, v2, v3, v4);
            }
        } else if (m_solidSimulationMethod == 3) {
            TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++) {
                const unsigned int v1 = tets[4 * i] + offset;
                const unsigned int v2 = tets[4 * i + 1] + offset;
                const unsigned int v3 = tets[4 * i + 2] + offset;
                const unsigned int v4 = tets[4 * i + 3] + offset;

                m_model.addStrainTetConstraint(v1, v2, v3, v4);
            }
        } else if (m_solidSimulationMethod == 4) {
            TetModel::ParticleMesh &mesh = m_model.getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++) {
                const unsigned int v[4] = {tets[4 * i] + offset, tets[4 * i + 1] + offset, tets[4 * i + 2] + offset,
                                           tets[4 * i + 3] + offset};
                // Important: Divide position correction by the number of clusters
                // which contain the vertex.
                const unsigned int nc[4] = {vTets[v[0]].m_numTets, vTets[v[1]].m_numTets, vTets[v[2]].m_numTets,
                                            vTets[v[3]].m_numTets};
                m_model.addShapeMatchingConstraint(4, v, nc);
            }
        }
    }
}
}  // namespace vox