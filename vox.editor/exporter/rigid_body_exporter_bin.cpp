//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.editor/exporter/rigid_body_exporter_bin.h"

#include "vox.base/file_system.h"
#include "vox.editor/scene_configuration.h"
#include "vox.sph/simulation.h"

using namespace vox;
using namespace vox::utility;

RigidBodyExporter_BIN::RigidBodyExporter_BIN(SimulatorBase* base) : ExporterBase(base) { m_isFirstFrame = true; }

RigidBodyExporter_BIN::~RigidBodyExporter_BIN() = default;

void RigidBodyExporter_BIN::init(const std::string& outputPath) {
    m_exportPath = FileSystem::normalizePath(outputPath + "/rigid_bodies");
}

void RigidBodyExporter_BIN::step(const unsigned int frame) {
    if (!m_active) return;

    std::string fileName = "rb_data_";
    fileName = fileName + std::to_string(frame) + ".bin";
    std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);

    writeRigidBodies(exportFileName);
}

void RigidBodyExporter_BIN::reset() { m_isFirstFrame = true; }

void RigidBodyExporter_BIN::setActive(const bool active) {
    ExporterBase::setActive(active);
    if (m_active) FileSystem::makeDirs(m_exportPath);
}

void RigidBodyExporter_BIN::writeRigidBodies(const std::string& fileName) {
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nBoundaryModels = sim->numberOfBoundaryModels();

    const std::string& sceneFile = SceneConfiguration::getCurrent()->getSceneFile();
    const utility::SceneLoader::Scene& scene = SceneConfiguration::getCurrent()->getScene();

    std::string scene_path = FileSystem::getFilePath(sceneFile);

    // check if we have a static model
    bool isStatic = true;
    for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) {
        BoundaryModel* bm = sim->getBoundaryModel(i);
        if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated()) {
            isStatic = false;
            break;
        }
    }

    BinaryFileWriter binWriter;
    if (m_isFirstFrame || !isStatic) binWriter.openFile(fileName);

    if (m_isFirstFrame) {
        binWriter.write(nBoundaryModels);

        for (auto boundaryModel : scene.boundaryModels) {
            std::string meshFileName = boundaryModel->meshFile;
            if (FileSystem::isRelativePath(meshFileName))
                meshFileName = FileSystem::normalizePath(scene_path + "/" + meshFileName);

            const std::string fileNameMesh = utility::FileSystem::getFileNameWithExt(meshFileName);
            binWriter.write(fileNameMesh);
            Eigen::Vector3f s = boundaryModel->scale.template cast<float>();
            binWriter.writeMatrix(s);
            std::string targetFilePath = m_exportPath + "/" + fileNameMesh;
            if (!utility::FileSystem::fileExists(targetFilePath)) {
                utility::FileSystem::copyFile(meshFileName, targetFilePath);
            }
            binWriter.write((char)boundaryModel->isWall);
            binWriter.writeMatrix(boundaryModel->color);
        }
    }

    if (m_isFirstFrame || !isStatic) {
        for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++) {
            BoundaryModel* bm = sim->getBoundaryModel(i);
            const Vector3r& x = bm->getRigidBodyObject()->getWorldSpacePosition();
            const Eigen::Vector3f x_f = x.template cast<float>();
            binWriter.writeMatrix(x_f);

            const Matrix3r& R = bm->getRigidBodyObject()->getWorldSpaceRotation();
            const Eigen::Matrix3f R_f = R.template cast<float>();
            // const Eigen::Matrix3f RT = R.transpose().template cast<float>();
            binWriter.writeMatrix(R_f);
        }
        binWriter.closeFile();
    }

    m_isFirstFrame = false;
}