//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.editor/position_based_dynamics_wrapper/pbd_boundary_simulator.h"

#include "vox.base/logging.h"
#include "vox.base/partio_reader_writer.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.editor/gui/simulator_gui_base.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_rigid_body.h"
#include "vox.editor/scene_configuration.h"
#include "vox.sph/emitter.h"
#include "vox.sph/sampling/surface_sampling.h"
#include "vox.sph/simulation.h"

using namespace std;
using namespace vox;
using namespace vox::utility;

PBDBoundarySimulator::PBDBoundarySimulator(SimulatorBase *base) {
    m_pbdWrapper = new PBDWrapper();
    m_base = base;
}

PBDBoundarySimulator::~PBDBoundarySimulator() { delete m_pbdWrapper; }

void PBDBoundarySimulator::init() {
    const utility::SceneLoader::Scene &scene = SceneConfiguration::getCurrent()->getScene();
    if (scene.sim2D) {
        LOGE("Dynamic boundaries are not supported in 2D simulations.")
        exit(1);
    }
}

void PBDBoundarySimulator::timeStep() {
    updateBoundaryForces();

    //////////////////////////////////////////////////////////////////////////
    // PBD
    //////////////////////////////////////////////////////////////////////////
    START_TIMING("SimStep - PBD")
    m_pbdWrapper->timeStep();
    STOP_TIMING_AVG

    Simulation *sim = Simulation::getCurrent();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
        m_base->updateBoundaryParticles(false);
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();
}

void PBDBoundarySimulator::initBoundaryData() {
    const std::string &sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
    const utility::SceneLoader::Scene &scene = SceneConfiguration::getCurrent()->getScene();

    // create additional rigid body information for emitters
    std::vector<PBDWrapper::RBData> additionalRBs;
    const std::string scene_path = FileSystem::getFilePath(sceneFile);
    // the last boundary models are the ones that were added for the emitters
    for (auto ed : scene.emitters) {
        PBDWrapper::RBData rb;
        rb.x = ed->x;
        rb.R = ed->rotation;
        rb.scale = Emitter::getSize(static_cast<Real>(ed->width), static_cast<Real>(ed->height), ed->type);
        rb.restitution = 0.6;
        rb.friction = 0.1;
        if (ed->type == 0) {
            rb.objFile = FileSystem::normalizePath(m_base->getExePath() + "/resources/emitter_boundary/EmitterBox.obj");
            rb.collisionType = 2;
        } else if (ed->type == 1) {
            rb.objFile =
                    FileSystem::normalizePath(m_base->getExePath() + "/resources/emitter_boundary/EmitterCylinder.obj");
            rb.collisionType = 5;
        }
        additionalRBs.push_back(rb);
    }
    m_pbdWrapper->readScene(sceneFile, additionalRBs);

    std::string scene_file_name = FileSystem::getFileName(sceneFile);
    const bool useCache = m_base->getUseParticleCaching();
    Simulation *sim = Simulation::getCurrent();

    string cachePath = scene_path + "/Cache";

    for (unsigned int i = 0; i < scene.boundaryModels.size(); i++) {
        string meshFileName = scene.boundaryModels[i]->meshFile;
        if (FileSystem::isRelativePath(meshFileName))
            meshFileName = FileSystem::normalizePath(scene_path + "/" + scene.boundaryModels[i]->meshFile);

        // check if mesh file has changed
        std::string md5FileName =
                FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(meshFileName) + ".md5");
        bool md5 = false;
        if (useCache) {
            string md5Str = FileSystem::getFileMD5(meshFileName);
            if (FileSystem::fileExists(md5FileName)) md5 = FileSystem::checkMD5(md5Str, md5FileName);
        }

        // if a samples file is given, use this one
        std::vector<Vector3r> boundaryParticles;

        SimulationModel &model = m_pbdWrapper->getSimulationModel();
        SimulationModel::RigidBodyVector &rigidBodies = model.getRigidBodies();
        auto *rb = new PBDRigidBody(rigidBodies[i]);
        RigidBodyGeometry &geo = rigidBodies[i]->getGeometry();
        utility::IndexedFaceMesh &mesh = geo.getMesh();
        VertexData &vd = geo.getVertexData();

        if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
            // if a samples file is given, use this one
            if (!scene.boundaryModels[i]->samplesFile.empty()) {
                string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
                PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(),
                                                  scene.boundaryModels[i]->scale[0], boundaryParticles);

                // transform back to local coordinates
                for (auto &boundaryParticle : boundaryParticles)
                    boundaryParticle = rb->getRotation().toRotationMatrix().transpose() *
                                       (rb->getWorldSpaceRotation() * boundaryParticle + rb->getWorldSpacePosition() -
                                        rb->getPosition());
            } else  // if no samples file is given, sample the surface model
            {
                // Cache sampling
                std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
                std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

                const string resStr = StringTools::real2String(scene.boundaryModels[i]->scale[0]) + "_" +
                                      StringTools::real2String(scene.boundaryModels[i]->scale[1]) + "_" +
                                      StringTools::real2String(scene.boundaryModels[i]->scale[2]);

                const string modeStr = "_m" + std::to_string(scene.boundaryModels[i]->samplingMode);
                const string particleFileName = FileSystem::normalizePath(
                        cachePath + "/" + mesh_file_name + "_db_" + StringTools::real2String(scene.particleRadius) +
                        "_" + resStr + modeStr + ".bgeo");

                // check MD5 if cache file is available
                bool foundCacheFile = false;

                if (useCache) foundCacheFile = FileSystem::fileExists(particleFileName);

                if (useCache && foundCacheFile && md5) {
                    PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(), 1.0,
                                                      boundaryParticles);
                    LOGI("Loaded cached boundary sampling: {}", particleFileName)
                }

                if (!useCache || !foundCacheFile || !md5) {
                    const auto samplePoissonDisk = [&]() {
                        LOGI("Poisson disk surface sampling of {}", meshFileName)
                        START_TIMING("Poisson disk sampling")
                        PoissonDiskSampling sampling;
                        sampling.sampleMesh(mesh.numVertices(), &vd.getPosition(0), mesh.numFaces(),
                                            mesh.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
                        STOP_TIMING_AVG
                    };
                    const auto sampleRegularTriangle = [&]() {
                        LOGI("Regular triangle surface sampling of {}", meshFileName)
                        START_TIMING("Regular triangle sampling")
                        RegularTriangleSampling sampling;
                        vox::RegularTriangleSampling::sampleMesh(mesh.numVertices(), &vd.getPosition(0),
                                                                 mesh.numFaces(), mesh.getFaces().data(),
                                                                 1.5f * scene.particleRadius, boundaryParticles);
                        STOP_TIMING_AVG
                    };
                    if (SurfaceSamplingMode::PoissonDisk == scene.boundaryModels[i]->samplingMode)
                        samplePoissonDisk();
                    else if (SurfaceSamplingMode::RegularTriangle == scene.boundaryModels[i]->samplingMode)
                        sampleRegularTriangle();
                    else {
                        LOGW("Unknown surface sampling method: {}", scene.boundaryModels[i]->samplingMode)
                        LOGW("Falling back to:")
                        sampleRegularTriangle();
                    }

                    // transform back to local coordinates
                    for (auto &boundaryParticle : boundaryParticles)
                        boundaryParticle = rb->getRotation().toRotationMatrix().transpose() *
                                           (boundaryParticle - rb->getPosition());

                    // Cache sampling
                    if (useCache && (FileSystem::makeDir(cachePath) == 0)) {
                        LOGI("Save particle sampling: {}", particleFileName)
                        PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(),
                                                           boundaryParticles.data(), nullptr, 0.0);
                    }
                }
            }

            auto *bm = new BoundaryModel_Akinci2012();
            bm->initModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
            sim->addBoundaryModel(bm);
        } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017) {
            auto *bm = new BoundaryModel_Koschier2017();
            bm->initModel(rb);
            sim->addBoundaryModel(bm);

            RigidBodyGeometry &geo = rigidBodies[i]->getGeometry();
            utility::IndexedFaceMesh &mesh = geo.getMesh();
            VertexData &vd = geo.getVertexData();
            // transform back
            std::vector<Vector3r> xLocal;
            xLocal.resize(vd.size());
            for (unsigned int j = 0; j < vd.size(); j++)
                xLocal[j] = rigidBodies[i]->getRotationMatrix().transpose() *
                            (vd.getPosition(j) - rigidBodies[i]->getPosition());

            m_base->initDensityMap(xLocal, mesh.getFaces(), scene.boundaryModels[i], md5, true, bm);
        } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019) {
            auto *bm = new BoundaryModel_Bender2019();
            bm->initModel(rb);
            sim->addBoundaryModel(bm);

            // transform back
            std::vector<Vector3r> xLocal;
            xLocal.resize(vd.size());
            for (unsigned int j = 0; j < vd.size(); j++)
                xLocal[j] = rb->getRotation().toRotationMatrix().transpose() * (vd.getPosition(j) - rb->getPosition());
            m_base->initVolumeMap(xLocal, mesh.getFaces(), scene.boundaryModels[i], md5, true, bm);
        }
        if (useCache && !md5) FileSystem::writeMD5File(meshFileName, md5FileName);
    }
}

void PBDBoundarySimulator::deferredInit() {
    Simulation *sim = Simulation::getCurrent();
    sim->performNeighborhoodSearchSort();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
        m_base->updateBoundaryParticles(true);
        Simulation::getCurrent()->updateBoundaryVolume();
    } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();

    m_pbdWrapper->initModel(TimeManager::getCurrent()->getTimeStepSize());
}

void PBDBoundarySimulator::reset() {
    //////////////////////////////////////////////////////////////////////////
    // PBD
    //////////////////////////////////////////////////////////////////////////
    m_pbdWrapper->reset();

    Simulation *sim = Simulation::getCurrent();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
        m_base->updateBoundaryParticles(true);
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();
}