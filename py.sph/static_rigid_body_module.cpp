//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "py.sph/common.h"
#include "vox.sph/static_rigid_body.h"

namespace py = pybind11;

void StaticRigidBodyModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Static Rigid Body
    // ---------------------------------------
    py::class_<vox::StaticRigidBody>(m_sub, "StaticRigidBody")
            .def(py::init<>([]() { return vox::StaticRigidBody(); }))
            .def("isDynamic", &vox::StaticRigidBody::isDynamic)
            .def("getMass", &vox::StaticRigidBody::getMass)
            .def("getPosition", &vox::StaticRigidBody::getPosition)
            .def("setPosition", &vox::StaticRigidBody::setPosition)
            .def("getWorldSpacePosition", &vox::StaticRigidBody::getWorldSpacePosition)
            .def("getVelocity", &vox::StaticRigidBody::getVelocity)
            .def("setVelocity", &vox::StaticRigidBody::setVelocity)
            .def("getRotation", &vox::StaticRigidBody::getRotation)
            .def("setRotation", &vox::StaticRigidBody::setRotation)
            .def("getWorldSpaceRotation", &vox::StaticRigidBody::getWorldSpaceRotation)
            .def("getAngularVelocity", &vox::StaticRigidBody::getAngularVelocity)
            .def("setAngularVelocity", &vox::StaticRigidBody::setAngularVelocity)
            .def("addForce", &vox::StaticRigidBody::addForce)
            .def("addTorque", &vox::StaticRigidBody::addTorque)

            .def("setWorldSpacePosition", &vox::StaticRigidBody::setWorldSpacePosition)
            .def("setWorldSpaceRotation", &vox::StaticRigidBody::setWorldSpaceRotation)
            .def("getGeometry", &vox::StaticRigidBody::getGeometry);
}
