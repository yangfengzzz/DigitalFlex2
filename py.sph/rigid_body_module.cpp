//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "py.sph/common.h"
#include "vox.sph/static_rigid_body.h"

namespace py = pybind11;

void RigidBodyModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Rigid Body Object
    // ---------------------------------------
    py::class_<vox::RigidBodyObject>(m_sub, "RigidBodyObject")
            .def("isDynamic", &vox::RigidBodyObject::isDynamic)
            .def("isAnimated", &vox::RigidBodyObject::isAnimated)
            .def("setIsAnimated", &vox::RigidBodyObject::setIsAnimated)
            .def("getMass", &vox::RigidBodyObject::getMass)
            .def("getPosition", &vox::RigidBodyObject::getPosition)
            .def("setPosition", &vox::RigidBodyObject::setPosition)
            .def("getWorldSpacePosition", &vox::RigidBodyObject::getWorldSpacePosition)
            .def("getVelocity", &vox::RigidBodyObject::getVelocity)
            .def("setVelocity", &vox::RigidBodyObject::setVelocity)
            .def("getRotation", &vox::RigidBodyObject::getRotation)
            .def("setRotation",
                 [](vox::RigidBodyObject& obj, const Vector4r& qVec) {
                     Quaternionr q;
                     q.coeffs() = qVec;
                     obj.setRotation(q);
                 })
            .def("getWorldSpaceRotation", &vox::RigidBodyObject::getWorldSpaceRotation)
            .def("getAngularVelocity", &vox::RigidBodyObject::getAngularVelocity)
            .def("setAngularVelocity", &vox::RigidBodyObject::setAngularVelocity)
            .def("addForce", &vox::RigidBodyObject::addForce)
            .def("addTorque", &vox::RigidBodyObject::addTorque)
            .def("getVertices", &vox::RigidBodyObject::getVertices)
            .def("getVertexNormals", &vox::RigidBodyObject::getVertexNormals)
            .def("getFaces", &vox::RigidBodyObject::getFaces)
            .def("updateMeshTransformation", &vox::RigidBodyObject::updateMeshTransformation);

    // ---------------------------------------
    // Class Static Rigid Body
    // ---------------------------------------
    py::class_<vox::StaticRigidBody, vox::RigidBodyObject>(m_sub, "StaticRigidBody")
            .def(py::init<>())
            .def("setWorldSpacePosition", &vox::StaticRigidBody::setWorldSpacePosition)
            .def("setWorldSpaceRotation", &vox::StaticRigidBody::setWorldSpaceRotation)
            .def("animate", &vox::StaticRigidBody::animate)
            .def("getGeometry", &vox::StaticRigidBody::getGeometry)
            .def("getVertexBuffer",
                 [](vox::StaticRigidBody& obj) -> py::memoryview {
                     auto vertices = obj.getVertices();
                     void* base_ptr = &vertices[0][0];
                     int num_vert = vertices.size();
                     return py::memoryview::from_buffer((Real*)base_ptr, {num_vert, 3},
                                                        {sizeof(Real) * 3, sizeof(Real)});
                 })
            .def("getFaceBuffer", [](vox::StaticRigidBody& obj) -> py::memoryview {
                auto faces = obj.getFaces();
                unsigned int* base_ptr = faces.data();
                return py::memoryview::from_buffer(base_ptr, {faces.size()}, {sizeof(unsigned int)});
            });
}
