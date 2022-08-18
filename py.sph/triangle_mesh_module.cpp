//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "common.h"
#include "vox.sph/simple_triangle_mesh.h"
#include "vox.sph/time_manager.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void TriangleMeshModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Time Manager
    // ---------------------------------------
    using Faces = vox::SimpleTriangleMesh::Faces;
    using Normals = vox::SimpleTriangleMesh::Normals;
    using Vertices = vox::SimpleTriangleMesh::Vertices;

    py::class_<vox::SimpleTriangleMesh>(m_sub, "TriangleMesh")
            .def(py::init<>())
            .def("release", &vox::SimpleTriangleMesh::release)
            .def("initMesh", &vox::SimpleTriangleMesh::initMesh)
            .def("addFace", overload_cast_<const unsigned int* const>()(&vox::SimpleTriangleMesh::addFace))
            .def("addFace", overload_cast_<const unsigned int* const>()(&vox::SimpleTriangleMesh::addFace))
            .def("addVertex", &vox::SimpleTriangleMesh::addVertex)

            .def("getFaces", (const Faces& (vox::SimpleTriangleMesh::*)() const)(&vox::SimpleTriangleMesh::getFaces))
            // .def("getFaces", (Faces & (vox::TriangleMesh::*)())(&vox::TriangleMesh::getFaces)) // TODO: wont work by
            // reference
            .def("getFaceNormals", (const Normals& (vox::SimpleTriangleMesh::*)() const)(&vox::SimpleTriangleMesh::getFaceNormals))
            // .def("getFaceNormals", (Normals & (vox::TriangleMesh::*)())(&vox::TriangleMesh::getFaceNormals)) // TODO:
            // wont work by reference
            .def("getVertexNormals",
                 (const Normals& (vox::SimpleTriangleMesh::*)() const)(&vox::SimpleTriangleMesh::getVertexNormals))
            // .def("getVertexNormals", (Normals & (vox::TriangleMesh::*)())(&vox::TriangleMesh::getVertexNormals)) //
            // TODO: wont work by reference
            .def("getVertices", (const Vertices& (vox::SimpleTriangleMesh::*)() const)(&vox::SimpleTriangleMesh::getVertices))
            .def("getVertexBuffer",
                 [](vox::SimpleTriangleMesh& obj) -> py::memoryview {
                     auto vertices = obj.getVertices();
                     void* base_ptr = &vertices[0][0];
                     int num_vert = obj.numVertices();
                     return py::memoryview::from_buffer((Real*)base_ptr, {num_vert, 3},
                                                        {sizeof(Real) * 3, sizeof(Real)});
                 })
            .def("getFaceBuffer",
                 [](vox::SimpleTriangleMesh& obj) -> py::memoryview {
                     auto faces = obj.getFaces();
                     void* base_ptr = faces.data();
                     int num_faces = obj.numFaces();
                     return py::memoryview::from_buffer((unsigned int*)base_ptr, {3 * num_faces},
                                                        {sizeof(unsigned int)});
                 })

            // .def("getVertices", (Vertices & (vox::TriangleMesh::*)())(&vox::TriangleMesh::getVertices)) // TODO: wont
            // work by reference

            .def("numVertices", &vox::SimpleTriangleMesh::numVertices)
            .def("numFaces", &vox::SimpleTriangleMesh::numFaces)

            .def("updateMeshTransformation", &vox::SimpleTriangleMesh::updateMeshTransformation)
            .def("updateNormals", &vox::SimpleTriangleMesh::updateNormals)
            .def("updateVertexNormals", &vox::SimpleTriangleMesh::updateVertexNormals);
}
