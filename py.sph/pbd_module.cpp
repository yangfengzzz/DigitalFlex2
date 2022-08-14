//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.editor/boundary_simulator.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_boundary_simulator.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_rigid_body.h"

namespace py = pybind11;

void PBDModule(py::module m) {
    auto m_sub = m.def_submodule("PBD");

    py::class_<vox::PBDWrapper>(m_sub, "PBDWrapper")
            .def(py::init<>())
            .def("reset", &vox::PBDWrapper::reset)
            .def("initModel", &vox::PBDWrapper::initModel)
            .def("readScene", &vox::PBDWrapper::readScene)
            .def("initTriangleModelConstraints", &vox::PBDWrapper::initTriangleModelConstraints)
            .def("initTetModelConstraints", &vox::PBDWrapper::initTetModelConstraints)
            .def("timeStep", &vox::PBDWrapper::timeStep)
            .def("updateVisModels", &vox::PBDWrapper::updateVisModels)
            .def("loadObj", &vox::PBDWrapper::loadObj)
            // .def("getSimulationModel", &PBDWrapper::getSimulationModel)  //TODO: make this work
            // .def("getCollisionDetection", &PBDWrapper::getCollisionDetection)  //TODO: make this work
            // .def("getTimeStepController", &PBDWrapper::getTimeStepController)  //TODO: make this work
            .def("getDampingCoeff", &vox::PBDWrapper::getDampingCoeff)
            .def("setDampingCoeff", &vox::PBDWrapper::setDampingCoeff)
            .def("getClothSimulationMethod", &vox::PBDWrapper::getClothSimulationMethod)
            .def("setClothSimulationMethod", &vox::PBDWrapper::setClothSimulationMethod)
            .def("getSolidSimulationMethod", &vox::PBDWrapper::getSolidSimulationMethod)
            .def("setSolidSimulationMethod", &vox::PBDWrapper::setSolidSimulationMethod)
            .def("getBendingMethod", &vox::PBDWrapper::getBendingMethod)
            .def("setBendingMethod", &vox::PBDWrapper::setBendingMethod);

    py::class_<vox::PBDBoundarySimulator, vox::BoundarySimulator>(m_sub, "PBDBoundarySimulator")
            .def(py::init<vox::SimulatorBase*>());
    //.def("getPBDWrapper", &SPH::PBDBoundarySimulator::getPBDWrapper, py::return_value_policy::reference_internal);
    ////TODO: make this work

    py::class_<vox::PBDRigidBody, vox::RigidBodyObject>(m_sub, "PBDRigidBody")
            .def(py::init<vox::RigidBody*>())
            .def("getVertices", &vox::PBDRigidBody::getVertices)
            .def("getFaces", &vox::PBDRigidBody::getFaces)
            .def("getVertexBuffer",
                 [](vox::PBDRigidBody& obj) -> py::memoryview {
                     const std::vector<Vector3r>& vertices = obj.getVertices();
                     void* base_ptr = const_cast<Real*>(&vertices[0][0]);
                     int num_vert = vertices.size();
                     return py::memoryview::from_buffer((Real*)base_ptr, {num_vert, 3},
                                                        {sizeof(Real) * 3, sizeof(Real)}, true);
                 })
            .def("getFaceBuffer", [](vox::PBDRigidBody& obj) -> py::memoryview {
                const std::vector<unsigned int>& faces = obj.getFaces();
                auto* base_ptr = const_cast<unsigned int*>(&faces[0]);
                return py::memoryview::from_buffer(base_ptr, {(int)faces.size()}, {sizeof(unsigned int)}, true);
            });
}
