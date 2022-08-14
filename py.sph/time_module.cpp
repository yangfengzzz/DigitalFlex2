//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/time_manager.h"
#include "vox.sph/time_step.h"

namespace py = pybind11;

void TimeModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Time Manager
    // ---------------------------------------
    py::class_<vox::TimeManager>(m_sub, "TimeManager")
            .def(py::init<>())
            .def_static("getCurrent", &vox::TimeManager::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &vox::TimeManager::setCurrent)
            .def_static("hasCurrent", &vox::TimeManager::hasCurrent)

            .def("getTime", &vox::TimeManager::getTime)
            .def("setTime", &vox::TimeManager::setTime)
            .def("getTimeStepSize", &vox::TimeManager::getTimeStepSize)
            .def("setTimeStepSize", &vox::TimeManager::setTimeStepSize)

            .def("saveState", &vox::TimeManager::saveState)
            .def("loadState", &vox::TimeManager::loadState);

    // ---------------------------------------
    // Abstract Class Time Step
    // ---------------------------------------
    py::class_<vox::TimeStep, vox::ParameterObject>(m_sub, "TimeStep")
            .def_readwrite_static("SOLVER_ITERATIONS", &vox::TimeStep::SOLVER_ITERATIONS)
            .def_readwrite_static("MIN_ITERATIONS", &vox::TimeStep::MIN_ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::TimeStep::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::TimeStep::MAX_ERROR)

            .def("step", &vox::TimeStep::step)
            .def("reset", &vox::TimeStep::reset)
            .def("init", &vox::TimeStep::init)
            .def("resize", &vox::TimeStep::resize)
            .def("emittedParticles", &vox::TimeStep::emittedParticles)
            .def("saveState", &vox::TimeStep::saveState)
            .def("loadState", &vox::TimeStep::loadState);
}