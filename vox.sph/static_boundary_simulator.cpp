//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/static_boundary_simulator.h"

#include "vox.base/logging.h"
#include "vox.sph/file_system.h"
#include "vox.sph/obj_loader.h"
#include "vox.sph/partio_reader_writer.h"
#include "vox.sph/scene_configuration.h"
#include "vox.sph/simulation.h"
#include "vox.sph/simulator_base.h"
#include "vox.sph/static_rigid_body.h"
#include "vox.sph/timing.h"
#include "vox.sph/utilities/surface_sampling.h"

using namespace std;
using namespace vox;
using namespace vox::utility;

StaticBoundarySimulator::StaticBoundarySimulator(SimulatorBase *base) { m_base = base; }

StaticBoundarySimulator::~StaticBoundarySimulator() = default;

void StaticBoundarySimulator::loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale) {
    std::vector<OBJLoader::Vec3f> x;
    std::vector<OBJLoader::Vec3f> normals;
    std::vector<MeshFaceIndices> faces;
    OBJLoader::Vec3f s = {(float)scale[0], (float)scale[1], (float)scale[2]};
    OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

    mesh.release();
    const auto nPoints = (unsigned int)x.size();
    const auto nFaces = (unsigned int)faces.size();
    mesh.initMesh(nPoints, nFaces);
    for (unsigned int i = 0; i < nPoints; i++) {
        mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
    }
    for (unsigned int i = 0; i < nFaces; i++) {
        // Reduce the indices by one
        int posIndices[3];
        for (int j = 0; j < 3; j++) {
            posIndices[j] = faces[i].posIndices[j] - 1;
        }

        mesh.addFace(&posIndices[0]);
    }

    LOGI("Number of triangles: {}", nFaces)
    LOGI("Number of vertices: {}", nPoints)
}

void StaticBoundarySimulator::initBoundaryData() {
    const std::string &sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
    const utility::SceneLoader::Scene &scene = SceneConfiguration::getCurrent()->getScene();
    std::string scene_path = FileSystem::getFilePath(sceneFile);
    std::string scene_file_name = FileSystem::getFileName(sceneFile);
    // no cache for 2D scenes
    // 2D sampling is fast, but storing it would require storing the transformation as well
    const bool useCache = m_base->getUseParticleCaching() && !scene.sim2D;
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

        auto *rb = new StaticRigidBody();
        rb->setIsAnimated(scene.boundaryModels[i]->isAnimated);
        TriangleMesh &geo = rb->getGeometry();
        loadObj(meshFileName, geo, scene.boundaryModels[i]->scale);

        std::vector<Vector3r> boundaryParticles;
        if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
            // if a samples file is given, use this one
            if (!scene.boundaryModels[i]->samplesFile.empty()) {
                string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
                PartioReaderWriter::readParticles(particleFileName, Vector3r::Zero(), Matrix3r::Identity(),
                                                  scene.boundaryModels[i]->scale[0], boundaryParticles);
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
                        cachePath + "/" + mesh_file_name + "_sb_" + StringTools::real2String(scene.particleRadius) +
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
                    if (!scene.sim2D) {
                        const auto samplePoissonDisk = [&]() {
                            LOGI("Poisson disk surface sampling of {}", meshFileName)
                            START_TIMING("Poisson disk sampling")
                            PoissonDiskSampling sampling;
                            sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(),
                                                geo.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
                            STOP_TIMING_AVG
                        };
                        const auto sampleRegularTriangle = [&]() {
                            LOGI("Regular triangle surface sampling of {}", meshFileName)
                            START_TIMING("Regular triangle sampling")
                            RegularTriangleSampling sampling;
                            vox::RegularTriangleSampling::sampleMesh(geo.numVertices(), geo.getVertices().data(),
                                                                     geo.numFaces(), geo.getFaces().data(),
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
                    } else {
                        LOGI("2D regular sampling of {}", meshFileName)
                        START_TIMING("2D regular sampling")
                        RegularSampling2D sampling;
                        vox::RegularSampling2D::sampleMesh(
                                Matrix3r::Identity(), Vector3r::Zero(), geo.numVertices(), geo.getVertices().data(),
                                geo.numFaces(), geo.getFaces().data(), 1.75f * scene.particleRadius, boundaryParticles);
                        STOP_TIMING_AVG
                    }

                    // Cache sampling
                    if (useCache && (FileSystem::makeDir(cachePath) == 0)) {
                        LOGI("Save particle sampling: {}", particleFileName)
                        PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(),
                                                           boundaryParticles.data(), nullptr, scene.particleRadius);
                    }
                }
            }
        }

        Quaternionr q(scene.boundaryModels[i]->rotation);
        rb->setPosition0(scene.boundaryModels[i]->translation);
        rb->setPosition(scene.boundaryModels[i]->translation);
        rb->setRotation0(q);
        rb->setRotation(q);

        if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
            auto *bm = new BoundaryModel_Akinci2012();
            bm->initModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
            sim->addBoundaryModel(bm);
        } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017) {
            auto *bm = new BoundaryModel_Koschier2017();
            bm->initModel(rb);
            sim->addBoundaryModel(bm);
            TriangleMesh &mesh = rb->getGeometry();
            m_base->initDensityMap(mesh.getVertices(), mesh.getFaces(), scene.boundaryModels[i], md5, false, bm);
        } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019) {
            auto *bm = new BoundaryModel_Bender2019();
            bm->initModel(rb);
            sim->addBoundaryModel(bm);
            TriangleMesh &mesh = rb->getGeometry();
            m_base->initVolumeMap(mesh.getVertices(), mesh.getFaces(), scene.boundaryModels[i], md5, false, bm);
        }
        if (useCache && !md5) FileSystem::writeMD5File(meshFileName, md5FileName);
        rb->updateMeshTransformation();
    }
}

void StaticBoundarySimulator::deferredInit() {
    Simulation *sim = Simulation::getCurrent();
    sim->performNeighborhoodSearchSort();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
        m_base->updateBoundaryParticles(true);
        Simulation::getCurrent()->updateBoundaryVolume();
    } else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();

#ifdef GPU_NEIGHBORHOOD_SEARCH
    // copy the particle data to the GPU
    sim->getNeighborhoodSearch()->update_point_sets();
#endif
}

void StaticBoundarySimulator::timeStep() {
    Simulation *sim = Simulation::getCurrent();
    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
        m_base->updateBoundaryParticles(false);
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();
}

void StaticBoundarySimulator::reset() {
    Simulation *sim = Simulation::getCurrent();
    for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) {
        BoundaryModel *bm = sim->getBoundaryModel(i);
        if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated()) {
            ((StaticRigidBody *)bm->getRigidBodyObject())->reset();
        }
    }

    if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
        m_base->updateBoundaryParticles(true);
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
        m_base->updateDMVelocity();
    else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
        m_base->updateVMVelocity();
}
