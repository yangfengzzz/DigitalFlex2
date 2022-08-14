//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.pbd/common.h"
#include "vox.base/logging.h"
#include "vox.pbd/obj_loader.h"
#include "vox.pbd/tet_gen_loader.h"
#include "vox.base/timing.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void UtilitiesModule(const py::module& m_sub) {
    using Faces = vox::utility::IndexedFaceMesh::Faces;
    using FaceNormals = vox::utility::IndexedFaceMesh::FaceNormals;
    using VertexNormals = vox::utility::IndexedFaceMesh::VertexNormals;
    using UVIndices = vox::utility::IndexedFaceMesh::UVIndices;
    using UVs = vox::utility::IndexedFaceMesh::UVs;

    py::class_<vox::utility::IndexedFaceMesh>(m_sub, "IndexedFaceMesh")
            .def(py::init<>())
            .def("release", &vox::utility::IndexedFaceMesh::release)
            .def("isClosed", &vox::utility::IndexedFaceMesh::isClosed)
            .def("initMesh", &vox::utility::IndexedFaceMesh::initMesh)
            .def("addFace", overload_cast_<const unsigned int* const>()(&vox::utility::IndexedFaceMesh::addFace))
            .def("addFace", overload_cast_<const int* const>()(&vox::utility::IndexedFaceMesh::addFace))
            .def("addUV", &vox::utility::IndexedFaceMesh::addUV)
            .def("addUVIndex", &vox::utility::IndexedFaceMesh::addUVIndex)

            .def("getFaces",
                 [](vox::utility::IndexedFaceMesh& mesh) -> py::memoryview {
                     const std::vector<unsigned int>& faces = mesh.getFaces();
                     auto* base_ptr = const_cast<unsigned int*>(&faces[0]);
                     return py::memoryview::from_buffer(base_ptr, {(int)faces.size()}, {sizeof(unsigned int)}, true);
                 })
            .def("getFaceNormals",
                 [](vox::utility::IndexedFaceMesh& mesh) -> py::memoryview {
                     const auto& n = mesh.getFaceNormals();
                     Real* base_ptr = const_cast<Real*>(&n[0][0]);
                     return py::memoryview::from_buffer(base_ptr, {(int)n.size(), 3}, {sizeof(Real) * 3, sizeof(Real)},
                                                        true);
                 })
            .def("getVertexNormals",
                 [](vox::utility::IndexedFaceMesh& mesh) -> py::memoryview {
                     const auto& n = mesh.getVertexNormals();
                     Real* base_ptr = const_cast<Real*>(&n[0][0]);
                     return py::memoryview::from_buffer(base_ptr, {(int)n.size(), 3}, {sizeof(Real) * 3, sizeof(Real)},
                                                        true);
                 })
            .def("getEdges", (const vox::utility::IndexedFaceMesh::Edges& (vox::utility::IndexedFaceMesh::*)()
                                      const)(&vox::utility::IndexedFaceMesh::getEdges))
            // .def("getEdges", (Edges &
            // (vox::utility::IndexedFaceMesh::*)())(&vox::utility::IndexedFaceMesh::getEdges)) //
            // TODO: wont work by reference
            .def("getFacesEdges", (const vox::utility::IndexedFaceMesh::FacesEdges& (vox::utility::IndexedFaceMesh::*)()
                                           const)(&vox::utility::IndexedFaceMesh::getFacesEdges))
            .def("getUVIndices", (const UVIndices& (vox::utility::IndexedFaceMesh::*)()
                                          const)(&vox::utility::IndexedFaceMesh::getUVIndices))
            .def("getUVs",
                 (const UVs& (vox::utility::IndexedFaceMesh::*)() const)(&vox::utility::IndexedFaceMesh::getUVs))
            .def("getVertexFaces",
                 (const vox::utility::IndexedFaceMesh::VerticesFaces& (vox::utility::IndexedFaceMesh::*)()
                          const)(&vox::utility::IndexedFaceMesh::getVertexFaces))
            .def("getVertexEdges",
                 (const vox::utility::IndexedFaceMesh::VerticesEdges& (vox::utility::IndexedFaceMesh::*)()
                          const)(&vox::utility::IndexedFaceMesh::getVertexEdges))
            .def("numVertices", &vox::utility::IndexedFaceMesh::numVertices)
            .def("numFaces", &vox::utility::IndexedFaceMesh::numFaces)
            .def("numEdges", &vox::utility::IndexedFaceMesh::numEdges)
            .def("numUVs", &vox::utility::IndexedFaceMesh::numUVs)
            .def("copyUVs", &vox::utility::IndexedFaceMesh::copyUVs)
            .def("getVerticesPerFace", &vox::utility::IndexedFaceMesh::getVerticesPerFace)
            .def("buildNeighbors", &vox::utility::IndexedFaceMesh::buildNeighbors)
            .def("updateNormalsVertexData", &vox::utility::IndexedFaceMesh::updateNormals<vox::VertexData>)
            .def("updateVertexNormalsVertexData", &vox::utility::IndexedFaceMesh::updateVertexNormals<vox::VertexData>)
            .def("updateNormalsParticleData", &vox::utility::IndexedFaceMesh::updateNormals<vox::ParticleData>)
            .def("updateVertexNormalsParticleData",
                 &vox::utility::IndexedFaceMesh::updateVertexNormals<vox::ParticleData>);

    py::class_<vox::utility::IndexedFaceMesh::Edge>(m_sub, "IndexedFaceMeshEdge")
            .def(py::init<>())
            .def_readwrite("m_face", &vox::utility::IndexedFaceMesh::Edge::m_face)
            .def_readwrite("m_vert", &vox::utility::IndexedFaceMesh::Edge::m_vert);

    py::class_<vox::utility::IndexedTetMesh>(m_sub, "IndexedTetMesh")
            .def(py::init<>())
            .def("release", &vox::utility::IndexedTetMesh::release)
            .def("initMesh", &vox::utility::IndexedTetMesh::initMesh)
            .def("addTet", overload_cast_<const unsigned int* const>()(&vox::utility::IndexedTetMesh::addTet))
            .def("addTet", overload_cast_<const int* const>()(&vox::utility::IndexedTetMesh::addTet))

            .def("getFaces",
                 [](vox::utility::IndexedTetMesh& mesh) -> py::memoryview {
                     const std::vector<unsigned int>& faces = mesh.getFaces();
                     auto* base_ptr = const_cast<unsigned int*>(&faces[0]);
                     return py::memoryview::from_buffer(base_ptr, {(int)faces.size()}, {sizeof(unsigned int)}, true);
                 })
            .def("getTets",
                 [](vox::utility::IndexedTetMesh& mesh) -> py::memoryview {
                     const std::vector<unsigned int>& tets = mesh.getTets();
                     auto* base_ptr = const_cast<unsigned int*>(&tets[0]);
                     return py::memoryview::from_buffer(base_ptr, {(int)tets.size()}, {sizeof(unsigned int)}, true);
                 })
            .def("getEdges", (const vox::utility::IndexedTetMesh::Edges& (vox::utility::IndexedTetMesh::*)()
                                      const)(&vox::utility::IndexedTetMesh::getEdges))
            .def("getFaceData", (const vox::utility::IndexedTetMesh::FaceData& (vox::utility::IndexedTetMesh::*)()
                                         const)(&vox::utility::IndexedTetMesh::getFaceData))
            .def("getTetData", (const vox::utility::IndexedTetMesh::TetData& (vox::utility::IndexedTetMesh::*)()
                                        const)(&vox::utility::IndexedTetMesh::getTetData))
            .def("getVertexTets", (const vox::utility::IndexedTetMesh::VerticesTets& (vox::utility::IndexedTetMesh::*)()
                                           const)(&vox::utility::IndexedTetMesh::getVertexTets))
            .def("getVertexFaces",
                 (const vox::utility::IndexedTetMesh::VerticesFaces& (vox::utility::IndexedTetMesh::*)()
                          const)(&vox::utility::IndexedTetMesh::getVertexFaces))
            .def("getVertexEdges",
                 (const vox::utility::IndexedTetMesh::VerticesEdges& (vox::utility::IndexedTetMesh::*)()
                          const)(&vox::utility::IndexedTetMesh::getVertexEdges))
            .def("numVertices", &vox::utility::IndexedTetMesh::numVertices)
            .def("numFaces", &vox::utility::IndexedTetMesh::numFaces)
            .def("numTets", &vox::utility::IndexedTetMesh::numTets)
            .def("numEdges", &vox::utility::IndexedTetMesh::numEdges)
            .def("buildNeighbors", &vox::utility::IndexedTetMesh::buildNeighbors);

    py::class_<vox::utility::IndexedTetMesh::Edge>(m_sub, "IndexedTetMeshEdge")
            .def(py::init<>())
            .def_readwrite("m_vert", &vox::utility::IndexedTetMesh::Edge::m_vert);

    py::class_<vox::utility::IndexedTetMesh::Face>(m_sub, "IndexedTetMeshFace")
            .def(py::init<>())
            .def_readwrite("m_tets", &vox::utility::IndexedTetMesh::Face::m_tets)
            .def_readwrite("m_edges", &vox::utility::IndexedTetMesh::Face::m_edges);

    py::class_<vox::utility::IndexedTetMesh::Tet>(m_sub, "IndexedTetMeshTet")
            .def(py::init<>())
            .def_readwrite("m_faces", &vox::utility::IndexedTetMesh::Tet::m_faces)
            .def_readwrite("m_edges", &vox::utility::IndexedTetMesh::Tet::m_edges);

    py::class_<vox::utility::MeshFaceIndices>(m_sub, "MeshFaceIndices")
            .def(py::init<>())
            .def_readwrite("posIndices", &vox::utility::MeshFaceIndices::posIndices)
            .def_readwrite("texIndices", &vox::utility::MeshFaceIndices::texIndices)
            .def_readwrite("normalIndices", &vox::utility::MeshFaceIndices::normalIndices);

    py::class_<vox::utility::OBJLoader>(m_sub, "OBJLoader")
            .def(py::init<>())
            .def("loadObj",
                 [](const std::string& filename, const vox::utility::OBJLoader::Vec3f& scale) {
                     std::vector<vox::utility::OBJLoader::Vec3f> x;
                     std::vector<vox::utility::OBJLoader::Vec3f> normals;
                     std::vector<vox::utility::OBJLoader::Vec2f> texCoords;
                     std::vector<vox::utility::MeshFaceIndices> faces;
                     vox::utility::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, scale);
                     return py::make_tuple(x, normals, texCoords, faces);
                 })
            .def("loadObjToMesh", [](const std::string& filename, const Vector3r& scale) {
                vox::VertexData vd;
                vox::utility::IndexedFaceMesh mesh;
                std::vector<vox::utility::OBJLoader::Vec3f> x;
                std::vector<vox::utility::OBJLoader::Vec3f> normals;
                std::vector<vox::utility::OBJLoader::Vec2f> texCoords;
                std::vector<vox::utility::MeshFaceIndices> faces;
                vox::utility::OBJLoader::Vec3f s = {(float)scale[0], (float)scale[1], (float)scale[2]};
                vox::utility::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

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
                return py::make_tuple(vd, mesh);
            });

    py::class_<vox::utility::TetGenLoader>(m_sub, "TetGenLoader")
            .def(py::init<>())
            .def("loadTetFile",
                 [](const std::string& filename) {
                     std::vector<Vector3r> x;
                     std::vector<unsigned int> tets;
                     vox::utility::TetGenLoader::loadTetFile(filename, x, tets);
                     return py::make_tuple(x, tets);
                 })
            .def("loadTetgenModel", [](const std::string& nodeFilename, const std::string& eleFilename) {
                std::vector<Vector3r> x;
                std::vector<unsigned int> tets;
                vox::utility::TetGenLoader::loadTetgenModel(nodeFilename, eleFilename, x, tets);
                return py::make_tuple(x, tets);
            });

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
            .def_static("stopTimingPrint",
                        []() {
                            static int timing_timerId = -1;
                            vox::utility::Timing::stopTiming(true);
                        })
            .def_static("stopTimingAvgPrint",
                        []() {
                            static int timing_timerId = -1;
                            vox::utility::Timing::stopTiming(true, timing_timerId);
                        })
            .def_static("stopTimingAvg",
                        []() {
                            static int timing_timerId = -1;
                            vox::utility::Timing::stopTiming(false, timing_timerId);
                        })
            .def_static("printAverageTimes", &vox::utility::Timing::printAverageTimes)
            .def_static("printTimeSums", &vox::utility::Timing::printTimeSums);
}
