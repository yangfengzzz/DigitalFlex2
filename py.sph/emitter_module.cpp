//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "py.sph/bind_pointer_vector.h"
#include "py.sph/common.h"
#include "vox.sph/emitter_system.h"

namespace py = pybind11;

void EmitterModule(const py::module &m_sub) {
    // ---------------------------------------
    // Emitter Class
    // ---------------------------------------
    py::class_<vox::Emitter>(m_sub, "Emitter")
            .def(py::init<>([](vox::FluidModel *model, const unsigned int width, const unsigned int height,
                               const Vector3r &pos, const Matrix3r &rotation, const Real velocity,
                               const unsigned int type = 0) {
                return vox::Emitter(model, width, height, pos, rotation, velocity, type);
            }))
            .def("emitParticles", &vox::Emitter::emitParticles)
            .def("emitParticlesCircle", &vox::Emitter::emitParticlesCircle)
            .def("getNextEmitTime", &vox::Emitter::getNextEmitTime)
            .def("setNextEmitTime", &vox::Emitter::setNextEmitTime)
            .def("setEmitStartTime", &vox::Emitter::setEmitStartTime)
            .def("setEmitEndTime", &vox::Emitter::setEmitEndTime)
            .def("getPosition", &vox::Emitter::getPosition)
            .def("setPosition", &vox::Emitter::setPosition)
            .def("getRotation", &vox::Emitter::getRotation)
            .def("setRotation", &vox::Emitter::setRotation)
            .def("getVelocity", &vox::Emitter::getVelocity)
            .def("setVelocity", &vox::Emitter::setVelocity)
            .def_static("getSize", &vox::Emitter::getSize)
            .def("step", &vox::Emitter::step)
            .def("reset", &vox::Emitter::reset)
            .def("saveState", &vox::Emitter::saveState)
            .def("loadState", &vox::Emitter::loadState);

    py::bind_pointer_vector<std::vector<vox::Emitter *>>(m_sub, "EmitterVector");

    // ---------------------------------------
    // Emitter System Class
    // ---------------------------------------
    py::class_<vox::EmitterSystem>(m_sub, "EmitterSystem")
            .def(py::init<>([](vox::FluidModel *model) { return vox::EmitterSystem(model); }))
            .def("enableReuseParticles", &vox::EmitterSystem::enableReuseParticles)
            .def("disableReuseParticles", &vox::EmitterSystem::disableReuseParticles)
            .def("addEmitter", &vox::EmitterSystem::addEmitter)
            .def("numEmitters", &vox::EmitterSystem::numEmitters)
            .def("getEmitters", &vox::EmitterSystem::getEmitters, py::return_value_policy::reference_internal)
            .def("numReusedParticles", &vox::EmitterSystem::numReusedParticles)
            .def("numEmittedParticles", &vox::EmitterSystem::numEmittedParticles)
            .def("step", &vox::EmitterSystem::step)
            .def("reset", &vox::EmitterSystem::reset)
            .def("saveState", &vox::EmitterSystem::saveState)
            .def("loadState", &vox::EmitterSystem::loadState);
}
