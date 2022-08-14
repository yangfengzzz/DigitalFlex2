//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include "py.pbd/common.h"
#include "vox.pbd/time_step.h"
#include "vox.pbd/time_step_controller.h"

namespace py = pybind11;

void TimeStepModule(const py::module& m_sub) {
    py::class_<vox::TimeStep, vox::ParameterObject>(m_sub, "TimeStep")
            //.def(py::init<>())
            .def("step", &vox::TimeStep::step)
            .def("reset", &vox::TimeStep::reset)
            .def("init", &vox::TimeStep::init)
            .def("setCollisionDetection", &vox::TimeStep::setCollisionDetection)
            .def("getCollisionDetection", &vox::TimeStep::getCollisionDetection,
                 py::return_value_policy::reference_internal);

    py::class_<vox::TimeStepController, vox::TimeStep>(m_sub, "TimeStepController")
            .def_readwrite_static("NUM_SUB_STEPS", &vox::TimeStepController::NUM_SUB_STEPS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::TimeStepController::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS_V", &vox::TimeStepController::MAX_ITERATIONS_V)
            .def_readwrite_static("VELOCITY_UPDATE_METHOD", &vox::TimeStepController::VELOCITY_UPDATE_METHOD)
            .def_readwrite_static("ENUM_VUPDATE_FIRST_ORDER", &vox::TimeStepController::ENUM_VUPDATE_FIRST_ORDER)
            .def_readwrite_static("ENUM_VUPDATE_SECOND_ORDER", &vox::TimeStepController::ENUM_VUPDATE_SECOND_ORDER)

            .def(py::init<>());
}