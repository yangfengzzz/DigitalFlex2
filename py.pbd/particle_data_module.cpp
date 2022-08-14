//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include "py.pbd/common.h"
#include "vox.pbd/simulation_model.h"

namespace py = pybind11;

void ParticleDataModule(const py::module& m_sub) {
    py::class_<vox::VertexData>(m_sub, "VertexData")
            .def(py::init<>())
            .def("addVertex", &vox::VertexData::addVertex)
            .def("getPosition",
                 (const Vector3r& (vox::VertexData::*)(const unsigned int) const)(&vox::VertexData::getPosition))
            .def("setPosition", &vox::VertexData::setPosition)
            .def("resize", &vox::VertexData::resize)
            .def("reserve", &vox::VertexData::reserve)
            .def("release", &vox::VertexData::release)
            .def("size", &vox::VertexData::size)
            //.def("getVertices", &vox::VertexData::getVertices, py::return_value_policy::reference);
            .def("getVertices", [](vox::VertexData& vd) -> py::memoryview {
                void* base_ptr = const_cast<Real*>(&(vd.getVertices())[0][0]);
                int num_vert = vd.size();
                return py::memoryview::from_buffer((Real*)base_ptr, {num_vert, 3}, {sizeof(Real) * 3, sizeof(Real)},
                                                   true);
            });

    py::class_<vox::ParticleData>(m_sub, "ParticleData")
            .def(py::init<>())
            .def("addVertex", &vox::ParticleData::addVertex)
            .def("getPosition",
                 (const Vector3r& (vox::ParticleData::*)(const unsigned int) const)(&vox::ParticleData::getPosition))
            .def("setPosition", &vox::ParticleData::setPosition)
            .def("getPosition0",
                 (const Vector3r& (vox::ParticleData::*)(const unsigned int) const)(&vox::ParticleData::getPosition0))
            .def("setPosition0", &vox::ParticleData::setPosition0)
            .def("getMass", (Real(vox::ParticleData::*)(const unsigned int) const)(&vox::ParticleData::getMass))
            .def("getInvMass", &vox::ParticleData::getInvMass)
            .def("setMass", &vox::ParticleData::setMass)
            .def("getVelocity",
                 (const Vector3r& (vox::ParticleData::*)(const unsigned int) const)(&vox::ParticleData::getVelocity))
            .def("setVelocity", &vox::ParticleData::setVelocity)
            .def("getAcceleration", (const Vector3r& (vox::ParticleData::*)(const unsigned int)
                                             const)(&vox::ParticleData::getAcceleration))
            .def("setAcceleration", &vox::ParticleData::setAcceleration)
            .def("getNumberOfParticles", &vox::ParticleData::getNumberOfParticles)
            .def("resize", &vox::ParticleData::resize)
            .def("reserve", &vox::ParticleData::reserve)
            .def("release", &vox::ParticleData::release)
            .def("size", &vox::ParticleData::size)
            //.def("getVertices", &vox::ParticleData::getVertices, py::return_value_policy::reference);
            .def("getVertices", [](vox::ParticleData& pd) -> py::memoryview {
                void* base_ptr = const_cast<Real*>(&(pd.getVertices())[0][0]);
                int num_vert = pd.getNumberOfParticles();
                return py::memoryview::from_buffer((Real*)base_ptr, {num_vert, 3}, {sizeof(Real) * 3, sizeof(Real)},
                                                   true);
            });
}