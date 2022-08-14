//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <array>
#include <iostream>
#include <vector>

#include "py.sph/bind_pointer_vector.h"
#include "py.sph/common.h"
#include "vox.editor/partio_reader_writer.h"
#include "vox.sph/binary_file_reader_writer.h"
#include "vox.sph/counting.h"
#include "vox.sph/file_system.h"
#include "vox.sph/obj_loader.h"
#include "vox.sph/utilities/gauss_quadrature.h"
#include "vox.sph/utilities/math_functions.h"
#include "vox.sph/utilities/matrix_free_solver.h"
#include "vox.sph/utilities/poisson_disk_sampling.h"
#include "vox.sph/utilities/regular_sampling_2d.h"
#include "vox.sph/utilities/scene_loader.h"
#include "vox.sph/utilities/sdf_functions.h"
#include "vox.sph/utilities/simple_quadrature.h"
#include "vox.sph/utilities/surface_sampling.h"
#include "vox.sph/utilities/volume_sampling.h"
#include "vox.sph/utilities/winding_numbers.h"

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

namespace py = pybind11;

using namespace pybind11::literals;

void UtilitiesModule(py::module m) {
    // Utilities submodule
    auto m_sub = m.def_submodule("Utilities");

    // ---------------------------------------
    // Binary File reader and writer
    // ---------------------------------------
    py::class_<vox::BinaryFileWriter>(m_sub, "BinaryFileWriter")
            .def(py::init<>())
            .def("openFile", &vox::BinaryFileWriter::openFile)
            .def("closeFile", &vox::BinaryFileWriter::closeFile);

    py::class_<vox::BinaryFileReader>(m_sub, "BinaryFileReader")
            .def(py::init<>())
            .def("openFile", &vox::BinaryFileReader::openFile)
            .def("closeFile", &vox::BinaryFileReader::closeFile);

    // ---------------------------------------
    // Counting class
    // ---------------------------------------
    py::class_<vox::utility::AverageCount>(m_sub, "AverageCount")
            .def(py::init<>())
            .def_readwrite("sum", &vox::utility::AverageCount::sum)
            .def_readwrite("numberOfCalls", &vox::utility::AverageCount::numberOfCalls);

    py::class_<vox::utility::Counting>(m_sub, "Counting")
            .def(py::init<>())
            .def_readwrite_static("m_averageCounts", &vox::utility::Counting::m_averageCounts)
            .def_static("reset", &vox::utility::Counting::reset)
            .def_static("increaseCounter", &vox::utility::Counting::increaseCounter)
            .def_static("printAverageCounts", &vox::utility::Counting::printAverageCounts)
            .def_static("printCounterSums", &vox::utility::Counting::printCounterSums);

    // TODO: check if init counting need to be implemented
    m_sub.def("INCREASE_COUNTER", &vox::utility::Counting::increaseCounter);

    // ---------------------------------------
    // File System utilities
    // ---------------------------------------
    py::class_<vox::utility::FileSystem>(m_sub, "FileSystem")
            .def(py::init<>())
            .def_static("getFilePath", vox::utility::FileSystem::getFilePath)
            .def_static("getFileName", vox::utility::FileSystem::getFileName)
            .def_static("getFileNameWithExt", vox::utility::FileSystem::getFileNameWithExt)
            .def_static("getFileExt", vox::utility::FileSystem::getFileExt)
            .def_static("isRelativePath", vox::utility::FileSystem::isRelativePath)
            .def_static("makeDir", vox::utility::FileSystem::makeDir)
            .def_static("makeDirs", vox::utility::FileSystem::makeDirs)
            .def_static("normalizePath", vox::utility::FileSystem::normalizePath)
            .def_static("fileExists", vox::utility::FileSystem::fileExists)
            .def_static("getProgramPath", vox::utility::FileSystem::getProgramPath)
            .def_static("copyFile", vox::utility::FileSystem::copyFile)
            .def_static("isFile", vox::utility::FileSystem::isFile)
            .def_static("isDirectory", vox::utility::FileSystem::isDirectory)
            .def_static("getFilesInDirectory", vox::utility::FileSystem::getFilesInDirectory)
            .def_static("getFileMD5", vox::utility::FileSystem::getFileMD5)
            .def_static("writeMD5File", vox::utility::FileSystem::writeMD5File)
            .def_static("checkMD5", vox::utility::FileSystem::checkMD5);

    // ---------------------------------------
    // Object loader necessary vector types
    // ---------------------------------------
    py::bind_vector<std::vector<std::array<float, 3>>>(m_sub, "VecVec3f");
    py::bind_vector<std::vector<std::array<float, 2>>>(m_sub, "VecVec2f");
    py::bind_vector<std::vector<vox::utility::MeshFaceIndices>>(m_sub, "VecMeshFaceIndices");

    // Todo: Bind missing attributes
    py::class_<vox::utility::MeshFaceIndices>(m_sub, "MeshFaceIndices")
            .def(py::init<>())
            .def("__repr__", [](const vox::utility::MeshFaceIndices &m) {
                std::stringstream s;
                s << "posIndices\n" << m.posIndices[0] << " " << m.posIndices[1] << " " << m.posIndices[2] << "\n";
                s << "texIndices\n" << m.texIndices[0] << " " << m.texIndices[1] << " " << m.texIndices[2] << "\n";
                s << "normalIndices\n"
                  << m.normalIndices[0] << " " << m.normalIndices[1] << " " << m.normalIndices[2] << "\n";
                return s.str();
            });

    py::class_<vox::utility::OBJLoader>(m_sub, "OBJLoader")
            .def(py::init<>())
            .def_static("loadObj", &vox::utility::OBJLoader::loadObj);

    // ---------------------------------------
    // Partio reader and writer
    // ---------------------------------------
    py::bind_vector<std::vector<Vector3r>>(m_sub, "VecVector3r",
                                           "std::vector<Vector3r> not to be confused with"
                                           "VecVec3f which ist std::vector<std::array<float,3>>");

    py::class_<vox::utility::PartioReaderWriter>(m_sub, "PartioReaderWriter")
            .def(py::init<>())
            .def_static("readParticle",
                        overload_cast_<const std::string &, const Vector3r &, const Matrix3r &, const Real,
                                       std::vector<Vector3r> &>()(&vox::utility::PartioReaderWriter::readParticles))
            .def_static("readParticle", overload_cast_<const std::string &, const Vector3r &, const Matrix3r &,
                                                       const Real, std::vector<Vector3r> &, std::vector<Vector3r> &>()(
                                                &vox::utility::PartioReaderWriter::readParticles))
            .def_static("readParticle",
                        overload_cast_<const std::string &, const Vector3r &, const Matrix3r &, const Real,
                                       std::vector<Vector3r> &, std::vector<Vector3r> &, Real &>()(
                                &vox::utility::PartioReaderWriter::readParticles))
            .def_static("writeParticle", &vox::utility::PartioReaderWriter::writeParticles);

    // ---------------------------------------
    // String tools
    // ---------------------------------------
    // TODO: String tools are omitted, because python is basically one big string tool

    // ---------------------------------------
    // System Info
    // ---------------------------------------
    // TODO: System info is also omitted. Will add if explicitly required by other classes

    // ---------------------------------------
    // Timing
    // ---------------------------------------
    // TODO: Timing and find a way for everything to be actually printed
    py::class_<vox::utility::TimingHelper>(m_sub, "TimingHelper")
            .def(py::init<>())
            .def_readwrite("start", &vox::utility::TimingHelper::start)
            .def_readwrite("name", &vox::utility::TimingHelper::name);

    py::class_<vox::utility::AverageTime>(m_sub, "AverageTime")
            .def(py::init<>())
            .def_readwrite("totalTime", &vox::utility::AverageTime::totalTime)
            .def_readwrite("counter", &vox::utility::AverageTime::counter)
            .def_readwrite("name", &vox::utility::AverageTime::name);

    py::class_<vox::utility::IDFactory>(m_sub, "IDFactory")
            .def(py::init<>())
            .def_static("getId", &vox::utility::IDFactory::getId);

    py::class_<vox::utility::Timing>(m_sub, "Timing")
            .def(py::init<>())
            .def_readwrite_static("m_dontPrintTimes", &vox::utility::Timing::m_dontPrintTimes)
            .def_readwrite_static("m_startCounter", &vox::utility::Timing::m_startCounter)
            .def_readwrite_static("m_stopCounter", &vox::utility::Timing::m_stopCounter)
            .def_readwrite_static("m_timingStack", &vox::utility::Timing::m_timingStack)
            .def_readwrite_static("m_averageTimes", &vox::utility::Timing::m_averageTimes)
            .def_static("reset", &vox::utility::Timing::reset)
            .def_static("startTiming", &vox::utility::Timing::startTiming)
            .def_static("stopTiming", overload_cast_<bool>()(&vox::utility::Timing::stopTiming))
            .def_static("stopTiming", overload_cast_<bool, int &>()(&vox::utility::Timing::stopTiming))
            .def_static("printAverageTimes", &vox::utility::Timing::printAverageTimes)
            .def_static("printTimeSums", &vox::utility::Timing::printTimeSums);

    // ---------------------------------------
    // Gauss Quadrature
    // ---------------------------------------
    py::class_<vox::GaussQuadrature>(m_sub, "GaussQuadrature")
            .def_static("integrate", &vox::GaussQuadrature::integrate)
            .def_static("exportSamples", &vox::GaussQuadrature::exportSamples);

    // ---------------------------------------
    // Math Functions
    // ---------------------------------------
    py::class_<vox::MathFunctions>(m_sub, "MathFunctions")
            .def_static("extractRotation", &vox::MathFunctions::extractRotation)
            .def_static("pseudoInverse", &vox::MathFunctions::pseudoInverse)
            .def_static("svdWithInversionHandling", &vox::MathFunctions::svdWithInversionHandling)
            .def_static("eigenDecomposition", &vox::MathFunctions::eigenDecomposition)
            .def_static("jacobiRotate", &vox::MathFunctions::jacobiRotate)
            .def_static("getOrthogonalVectors", &vox::MathFunctions::getOrthogonalVectors);

    // ---------------------------------------
    // Matrix Replacement TODO: Pretty sure the user doesnt need the matrix replacement solved in the python API
    // ---------------------------------------
    py::class_<vox::MatrixReplacement>(m_sub, "MatrixReplacement");

    // ---------------------------------------
    // Poisson Disk Sampling
    // ---------------------------------------
    py::class_<vox::PoissonDiskSampling>(m_sub, "PoissonDiskSampling")
            .def(py::init<>())
            .def_static("floor", &vox::PoissonDiskSampling::floor)
            .def("sampleMesh", &vox::PoissonDiskSampling::sampleMesh);

    using CellPos = Eigen::Matrix<int, 3, 1, Eigen::DontAlign>;
    py::class_<vox::PoissonDiskSampling::InitialPointInfo>(m_sub, "InitialPointInfo")
            .def(py::init<>())
            .def(py::init<CellPos, Vector3r, unsigned int>())
            .def_readwrite("cP", &vox::PoissonDiskSampling::InitialPointInfo::cP)
            .def_readwrite("pos", &vox::PoissonDiskSampling::InitialPointInfo::pos)
            .def_readwrite("ID", &vox::PoissonDiskSampling::InitialPointInfo::ID);

    py::class_<vox::PoissonDiskSampling::HashEntry>(m_sub, "HashEntry")
            .def(py::init<>())
            .def_readwrite("samples", &vox::PoissonDiskSampling::HashEntry::samples)
            .def_readwrite("startIndex", &vox::PoissonDiskSampling::HashEntry::startIndex);

    // ---------------------------------------
    // Regular Sampling 2D
    // ---------------------------------------
    py::class_<vox::RegularSampling2D>(m_sub, "RegularSampling2D")
            .def_static("sampleMesh", &vox::RegularSampling2D::sampleMesh);

    // ---------------------------------------
    // Scene Loader
    // ---------------------------------------
    py::class_<vox::utility::SceneLoader>(m_sub, "SceneLoader")
            .def(py::init<>())
            .def("readScene", &vox::utility::SceneLoader::readScene);

    auto m_sub_sub = m_sub.def_submodule("SceneLoaderStructs");
    using SceneInfo = vox::utility::SceneLoader;
    py::class_<SceneInfo::Box>(m_sub_sub, "Box")
            .def(py::init<Vector3r, Vector3r>(), "minX"_a = Vector3r::Zero(), "maxX"_a = Vector3r::Ones())
            .def_readwrite("m_minX", &SceneInfo::Box::m_minX)
            .def_readwrite("m_maxX", &SceneInfo::Box::m_maxX);

    py::class_<SceneInfo::BoundaryData>(m_sub_sub, "BoundaryData")
            // .def(py::init<>())
            .def(py::init<std::string, std::string, Vector3r, Matrix3r, Vector3r, Real, bool, bool,
                          Eigen::Matrix<float, 4, 1, Eigen::DontAlign>, void *, std::string, bool, Real,
                          Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>, unsigned int, bool>(),
                 "samplesFile"_a = "", "meshFile"_a = "", "translation"_a = Vector3r::Zero(),
                 "rotation"_a = Matrix3r::Identity(), "scale"_a = Vector3r::Ones(), "density"_a = 1000.,
                 "isDynamic"_a = false, "isWall"_a = false, "color"_a = Eigen::Vector4f(1.f, 0.f, 0.f, 0.f),
                 "rigidBody"_a = nullptr, "mapFile"_a = "", "mapInvert"_a = false, "mapThickness"_a = 0.0,
                 "mapResolution"_a = Eigen::Matrix<unsigned int, 3, 1>(20, 20, 20), "samplingMode"_a = 0,
                 "isAnimated"_a = false)
            .def_readwrite("samplesFile", &SceneInfo::BoundaryData::samplesFile)
            .def_readwrite("meshFile", &SceneInfo::BoundaryData::meshFile)
            .def_readwrite("translation", &SceneInfo::BoundaryData::translation)
            .def_readwrite("rotation", &SceneInfo::BoundaryData::rotation)
            .def_readwrite("scale", &SceneInfo::BoundaryData::scale)
            .def_readwrite("density", &SceneInfo::BoundaryData::density)
            .def_readwrite("dynamic", &SceneInfo::BoundaryData::dynamic)
            .def_readwrite("isWall", &SceneInfo::BoundaryData::isWall)
            .def_readwrite("color", &SceneInfo::BoundaryData::color)
            .def_readwrite("rigidBody",
                           &SceneInfo::BoundaryData::rigidBody)  // TODO: find out if void pointer is problem
            .def_readwrite("mapFile", &SceneInfo::BoundaryData::mapFile)
            .def_readwrite("mapInvert", &SceneInfo::BoundaryData::mapInvert)
            .def_readwrite("mapThickness", &SceneInfo::BoundaryData::mapThickness)
            .def_readwrite("mapResolution", &SceneInfo::BoundaryData::mapResolution)
            .def_readwrite("samplingMode", &SceneInfo::BoundaryData::samplingMode)
            .def_readwrite("isAnimated", &SceneInfo::BoundaryData::isAnimated);

    py::class_<SceneInfo::FluidData>(m_sub_sub, "FluidData")
            .def(py::init<std::string, std::string, std::string, Vector3r, Matrix3r, Vector3r, Vector3r, Vector3r,
                          unsigned char, bool, std::array<unsigned int, 3>>(),
                 "id"_a = "Fluid", "visMeshFile"_a = "", "samplesFile"_a,
                 "translation"_a =
                         Vector3r::Zero(),  // TODO: samples file here has no default because it needs to be provided
                 "rotation"_a = Matrix3r::Identity(), "scale"_a = Vector3r::Ones(),
                 "initialVelocity"_a = Vector3r::Zero(), "initialAngularVelocity"_a = Vector3r::Zero(), "mode"_a = 0,
                 "invert"_a = false, "resolutionSDF"_a = std::array<unsigned int, 3>({20, 20, 20}))
            .def_readwrite("id", &SceneInfo::FluidData::id)
            .def_readwrite("visMeshFile", &SceneInfo::FluidData::visMeshFile)
            .def_readwrite("samplesFile", &SceneInfo::FluidData::samplesFile)
            .def_readwrite("translation", &SceneInfo::FluidData::translation)
            .def_readwrite("rotation", &SceneInfo::FluidData::rotation)
            .def_readwrite("scale", &SceneInfo::FluidData::scale)
            .def_readwrite("initialVelocity", &SceneInfo::FluidData::initialVelocity)
            .def_readwrite("initialAngularVelocity", &SceneInfo::FluidData::initialAngularVelocity)
            .def_readwrite("mode", &SceneInfo::FluidData::mode)
            .def_readwrite("invert", &SceneInfo::FluidData::invert)
            .def_readwrite("resolutionSDF", &SceneInfo::FluidData::resolutionSDF);

    py::class_<SceneInfo::FluidBlock>(m_sub_sub, "FluidBlock")
            .def(py::init<std::string, std::string, SceneInfo::Box, unsigned char, Vector3r, Vector3r>(),
                 "id"_a = "Fluid", "visMeshFile"_a = "", "box"_a = SceneInfo::Box({Vector3r::Zero(), Vector3r::Ones()}),
                 "mode"_a = 0, "initialVelocity"_a = Vector3r::Zero(), "initialVelocity"_a = Vector3r::Zero())
            .def_readwrite("id", &SceneInfo::FluidBlock::id)
            .def_readwrite("visMeshFile", &SceneInfo::FluidBlock::visMeshFile)
            .def_readwrite("box", &SceneInfo::FluidBlock::box)
            .def_readwrite("mode", &SceneInfo::FluidBlock::mode)
            .def_readwrite("initialVelocity", &SceneInfo::FluidBlock::initialVelocity)
            .def_readwrite("initialAngularVelocity", &SceneInfo::FluidBlock::initialAngularVelocity);

    py::class_<SceneInfo::EmitterData>(m_sub_sub, "EmitterData")
            .def(py::init<std::string, unsigned int, unsigned int, Vector3r, Real, Matrix3r, Real, Real,
                          unsigned int>(),
                 "id"_a = "Fluid", "width"_a = 5, "height"_a = 5, "x"_a = Vector3r::Zero(), "velocity"_a = 1,
                 "rotation"_a = Matrix3r::Identity(), "emitStartTime"_a = 0,
                 "emitEndTime"_a = std::numeric_limits<Real>::max(), "type"_a = 0)
            .def_readwrite("id", &SceneInfo::EmitterData::id)
            .def_readwrite("width", &SceneInfo::EmitterData::width)
            .def_readwrite("height", &SceneInfo::EmitterData::height)
            .def_readwrite("x", &SceneInfo::EmitterData::x)
            .def_readwrite("velocity", &SceneInfo::EmitterData::velocity)
            .def_readwrite("rotation", &SceneInfo::EmitterData::rotation)
            .def_readwrite("emitStartTime", &SceneInfo::EmitterData::emitStartTime)
            .def_readwrite("emitEndTime", &SceneInfo::EmitterData::emitEndTime)
            .def_readwrite("type", &SceneInfo::EmitterData::type);

    py::class_<SceneInfo::AnimationFieldData>(m_sub_sub, "AnimationFieldData")
            .def(py::init<std::string, std::string, std::string, std::string, unsigned int, Vector3r, Matrix3r,
                          Vector3r, Real, Real>(),
                 "particleFieldName"_a = "", "expressionX"_a = "", "expressionY"_a = "", "expressionZ"_a = "",
                 "shapeType"_a = 0, "x"_a = Vector3r::Zero(), "rotation"_a = Matrix3r::Identity(),
                 "scale"_a = Vector3r::Ones(), "startTime"_a = 0, "endTime"_a = std::numeric_limits<Real>::max())
            .def_readwrite("particleFieldName", &SceneInfo::AnimationFieldData::particleFieldName)
            // .def_readwrite("expression", &SceneInfo::AnimationFieldData::expression) // TODO: bind this type manually
            .def_readwrite("shapeType", &SceneInfo::AnimationFieldData::shapeType)
            .def_readwrite("x", &SceneInfo::AnimationFieldData::x)
            .def_readwrite("rotation", &SceneInfo::AnimationFieldData::rotation)
            .def_readwrite("scale", &SceneInfo::AnimationFieldData::scale)
            .def_readwrite("startTime", &SceneInfo::AnimationFieldData::startTime)
            .def_readwrite("endTime", &SceneInfo::AnimationFieldData::endTime);

    py::class_<SceneInfo::MaterialData>(m_sub_sub, "MaterialData")
            .def(py::init<>())
            .def(py::init<std::string, std::string, unsigned int, Real, Real, unsigned int, bool, Vector3r, Vector3r>(),
                 "id"_a, "colorField"_a = "velocity", "colorMapType"_a = 1, "minVal"_a = 0.0,
                 "maxVal"_a = 10.0,  // TODO: an id has to be provided
                 "maxEmitterParticles"_a = 10000, "emitterReuseParticles"_a = false,
                 "emitterBoxMin"_a = Vector3r(-1.0, -1.0, -1.0), "emitterBoxMax"_a = Vector3r(1.0, 1.0, 1.0))
            .def_readwrite("id", &SceneInfo::MaterialData::id)
            .def_readwrite("colorField", &SceneInfo::MaterialData::colorField)
            .def_readwrite("colorMapType", &SceneInfo::MaterialData::colorMapType)
            .def_readwrite("minVal", &SceneInfo::MaterialData::minVal)
            .def_readwrite("maxVal", &SceneInfo::MaterialData::maxVal)
            .def_readwrite("maxEmitterParticles", &SceneInfo::MaterialData::maxEmitterParticles)
            .def_readwrite("emitterReuseParticles", &SceneInfo::MaterialData::emitterReuseParticles)
            .def_readwrite("emitterBoxMin", &SceneInfo::MaterialData::emitterBoxMin)
            .def_readwrite("emitterBoxMax", &SceneInfo::MaterialData::emitterBoxMax);

    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::BoundaryData *>>(m_sub, "BoundaryDataVector");
    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::FluidData *>>(m_sub, "FluidDataVector");
    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::FluidBlock *>>(m_sub, "FluidBlockVector");
    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::EmitterData *>>(m_sub, "EmitterDataVector");
    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::AnimationFieldData *>>(m_sub,
                                                                                          "AnimationFieldDataVector");
    py::bind_pointer_vector<std::vector<vox::utility::SceneLoader::MaterialData *>>(m_sub, "MaterialData");

    py::class_<SceneInfo::Scene>(m_sub_sub, "Scene")
            .def(py::init<>())
            .def_readwrite("boundaryModels", &SceneInfo::Scene::boundaryModels)
            .def_readwrite("fluidModels", &SceneInfo::Scene::fluidModels)
            .def_readwrite("fluidBlocks", &SceneInfo::Scene::fluidBlocks)
            .def_readwrite("emitters", &SceneInfo::Scene::emitters)
            .def_readwrite("animatedFields", &SceneInfo::Scene::animatedFields)
            .def_readwrite("materials", &SceneInfo::Scene::materials)
            .def_readwrite("particleRadius", &SceneInfo::Scene::particleRadius)
            .def_readwrite("sim2D", &SceneInfo::Scene::sim2D)
            .def_readwrite("timeStepSize", &SceneInfo::Scene::timeStepSize)
            .def_readwrite("camPosition", &SceneInfo::Scene::camPosition)
            .def_readwrite("camLookat", &SceneInfo::Scene::camLookat);

    // ---------------------------------------
    // SDF Functions TODO: implement discregrid
    // ---------------------------------------
    py::class_<vox::CubicLagrangeDiscreteGrid>(m_sub, "DiscreteGrid");

    py::class_<vox::utility::SDFFunctions>(m_sub, "SDFFunctions")
            .def_static("generateSDF", &vox::utility::SDFFunctions::generateSDF)
            .def_static("computeBoundingBox", &vox::utility::SDFFunctions::computeBoundingBox)
            .def_static("distance", overload_cast_<vox::CubicLagrangeDiscreteGrid *, const Vector3r &, const Real,
                                                   Vector3r &, Vector3r &>()(&vox::utility::SDFFunctions::distance))
            .def_static("distance", overload_cast_<vox::CubicLagrangeDiscreteGrid *, const Vector3r &, const Real>()(
                                            &vox::utility::SDFFunctions::distance));

    // ---------------------------------------
    // Simple Quadrature
    // ---------------------------------------
    py::class_<vox::SimpleQuadrature>(m_sub, "SimpleQuadrature")
            .def_readwrite_static("m_samplePoints", &vox::SimpleQuadrature::m_samplePoints)
            .def_readwrite_static("m_volume", &vox::SimpleQuadrature::m_volume)

            .def_static("determineSamplePointsInSphere", &vox::SimpleQuadrature::determineSamplePointsInSphere)
            .def_static("determineSamplePointsInCircle", &vox::SimpleQuadrature::determineSamplePointsInCircle)
            .def_static("integrate", &vox::SimpleQuadrature::integrate);

    // ---------------------------------------
    // Surface Sampling Modes
    // ---------------------------------------
    py::enum_<vox::SurfaceSamplingMode>(m_sub, "SurfaceSamplingMode")
            .value("PoissonDisk", vox::SurfaceSamplingMode::PoissonDisk)
            .value("RegularTriangle", vox::SurfaceSamplingMode::RegularTriangle)
            .value("Regular2D", vox::SurfaceSamplingMode::Regular2D);

    // ---------------------------------------
    // Volume Sampling
    // ---------------------------------------
    py::class_<vox::utility::VolumeSampling>(m_sub, "VolumeSampling")
            .def_static("sampleMesh", vox::utility::VolumeSampling::sampleMesh);

    // ---------------------------------------
    // Winding Numbers
    // ---------------------------------------
    py::class_<vox::utility::WindingNumbers>(m_sub, "WindingNumbers")
            .def_static("computeGeneralizedWindingNumber",
                        overload_cast_<const Vector3r &, const Vector3r &, const Vector3r &, const Vector3r &>()(
                                &vox::utility::WindingNumbers::computeGeneralizedWindingNumber))
            .def_static("computeGeneralizedWindingNumber",
                        overload_cast_<const Vector3r &, const vox::TriangleMesh &>()(
                                &vox::utility::WindingNumbers::computeGeneralizedWindingNumber));
}
