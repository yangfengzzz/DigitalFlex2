//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>

#include "py.pbd/common.h"

namespace py = pybind11;
using namespace vox;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

#define GET_QUAT_FCT_CONVERT(name) .def(#name, [](vox::RigidBody& obj) { return obj.name().coeffs(); })
#define SET_QUAT_FCT_CONVERT(name)                              \
    .def(#name, [](vox::RigidBody& obj, const Vector4r& qVec) { \
        Quaternionr q;                                          \
        q.coeffs() = qVec;                                      \
        obj.name(q);                                            \
    })

void RigidBodyModule(const py::module& m_sub) {
    py::class_<vox::RigidBodyGeometry>(m_sub, "RigidBodyGeometry")
            .def(py::init<>())
            .def("getMesh", &vox::RigidBodyGeometry::getMesh)
            .def("getVertexData",
                 (const vox::VertexData& (vox::RigidBodyGeometry::*)() const)(&vox::RigidBodyGeometry::getVertexData))
            .def("getVertexDataLocal", (const vox::VertexData& (vox::RigidBodyGeometry::*)()
                                                const)(&vox::RigidBodyGeometry::getVertexDataLocal))
            .def("initMesh", &vox::RigidBodyGeometry::initMesh)
            .def("updateMeshTransformation", &vox::RigidBodyGeometry::updateMeshTransformation)
            .def("updateMeshNormals", &vox::RigidBodyGeometry::updateMeshNormals);

    py::class_<vox::RigidBody>(m_sub, "RigidBody")
            .def(py::init<>())
            .def("initBody",
                 [](vox::RigidBody& obj, const Real density, const Vector3r& x, const Vector4r& qVec,
                    const vox::VertexData& vertices, const vox::utility::IndexedFaceMesh& mesh, const Vector3r& scale) {
                     Quaternionr q;
                     q.coeffs() = qVec;
                     obj.initBody(density, x, q, vertices, mesh, scale);
                 })
            .def("initBody",
                 [](vox::RigidBody& obj, const Real mass, const Vector3r& x, const Vector3r& inertiaTensor,
                    const Vector4r& qVec, const vox::VertexData& vertices, const vox::utility::IndexedFaceMesh& mesh,
                    const Vector3r& scale) {
                     Quaternionr q;
                     q.coeffs() = qVec;
                     obj.initBody(mass, x, inertiaTensor, q, vertices, mesh, scale);
                 })
            .def("reset", &vox::RigidBody::reset)
            .def("updateInverseTransformation", &vox::RigidBody::updateInverseTransformation)
            .def("rotationUpdated", &vox::RigidBody::rotationUpdated)
            .def("updateInertiaW", &vox::RigidBody::updateInertiaW)
            .def("determineMassProperties", &vox::RigidBody::determineMassProperties)
            .def("getTransformationR", &vox::RigidBody::getTransformationR)
            .def("getTransformationV1", &vox::RigidBody::getTransformationV1)
            .def("getTransformationV2", &vox::RigidBody::getTransformationV2)
            .def("getTransformationRXV1", &vox::RigidBody::getTransformationRXV1)
            .def("getMass", (const Real& (vox::RigidBody::*)() const)(&vox::RigidBody::getMass))
            .def("setMass", &vox::RigidBody::setMass)
            .def("getInvMass", &vox::RigidBody::getInvMass)
            .def("getPosition", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getPosition))
            .def("setPosition", &vox::RigidBody::setPosition)
            .def("getLastPosition", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getLastPosition))
            .def("setLastPosition", &vox::RigidBody::setLastPosition)
            .def("getOldPosition", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getOldPosition))
            .def("setOldPosition", &vox::RigidBody::setOldPosition)
            .def("getPosition0", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getPosition0))
            .def("setPosition0", &vox::RigidBody::setPosition0)
            .def("getPositionInitial_MAT",
                 (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getPositionInitial_MAT))
            .def("setPositionInitial_MAT", &vox::RigidBody::setPositionInitial_MAT)
            .def("getVelocity", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getVelocity))
            .def("setVelocity", &vox::RigidBody::setVelocity)
            .def("getVelocity0", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getVelocity0))
            .def("setVelocity0", &vox::RigidBody::setVelocity0)
            .def("getAcceleration", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getAcceleration))
            .def("setAcceleration", &vox::RigidBody::setAcceleration)
            .def("getInertiaTensor", &vox::RigidBody::getInertiaTensor)
            .def("setInertiaTensor", &vox::RigidBody::setInertiaTensor)
            .def("getInertiaTensorW", (const Matrix3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getInertiaTensorW))
            .def("getInertiaTensorInverse", &vox::RigidBody::getInertiaTensorInverse)
            .def("getInertiaTensorInverseW",
                 (const Matrix3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getInertiaTensorInverseW))
            .def("setInertiaTensorInverseW", &vox::RigidBody::setInertiaTensorInverseW)
                    GET_QUAT_FCT_CONVERT(getRotation) SET_QUAT_FCT_CONVERT(setRotation)
                            GET_QUAT_FCT_CONVERT(getLastRotation) SET_QUAT_FCT_CONVERT(setLastRotation)
                                    GET_QUAT_FCT_CONVERT(getOldRotation) SET_QUAT_FCT_CONVERT(setOldRotation)
                                            GET_QUAT_FCT_CONVERT(getRotation0) SET_QUAT_FCT_CONVERT(setRotation0)
                                                    GET_QUAT_FCT_CONVERT(getRotationMAT)
                                                            SET_QUAT_FCT_CONVERT(setRotationMAT)
                                                                    GET_QUAT_FCT_CONVERT(getRotationInitial)
                                                                            SET_QUAT_FCT_CONVERT(setRotationInitial)
            .def("getRotationMatrix", (const Matrix3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getRotationMatrix))
            .def("setRotationMatrix", &vox::RigidBody::setRotationMatrix)
            .def("getAngularVelocity",
                 (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getAngularVelocity))
            .def("setAngularVelocity", &vox::RigidBody::setAngularVelocity)
            .def("getAngularVelocity0",
                 (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getAngularVelocity0))
            .def("setAngularVelocity0", &vox::RigidBody::setAngularVelocity0)
            .def("getTorque", (const Vector3r& (vox::RigidBody::*)() const)(&vox::RigidBody::getTorque))
            .def("setTorque", &vox::RigidBody::setTorque)
            .def("getRestitutionCoeff", &vox::RigidBody::getRestitutionCoeff)
            .def("setRestitutionCoeff", &vox::RigidBody::setRestitutionCoeff)
            .def("getFrictionCoeff", &vox::RigidBody::getFrictionCoeff)
            .def("setFrictionCoeff", &vox::RigidBody::setFrictionCoeff)
            .def("getGeometry", &vox::RigidBody::getGeometry);
}