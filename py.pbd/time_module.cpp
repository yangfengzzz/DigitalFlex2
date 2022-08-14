//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.pbd/common.h"
#include "vox.pbd/time_manager.h"

namespace py = pybind11;

namespace vox {

void TimeModule(py::module& m_sub) {
    // ---------------------------------------
    // Class Time Manager
    // ---------------------------------------
    py::class_<TimeManager>(m_sub, "TimeManager")
            .def(py::init<>())
            .def_static("getCurrent", &TimeManager::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &TimeManager::setCurrent)
            .def_static("hasCurrent", &TimeManager::hasCurrent)

            .def("getTime", &TimeManager::getTime)
            .def("setTime", &TimeManager::setTime)
            .def("getTimeStepSize", &TimeManager::getTimeStepSize)
            .def("setTimeStepSize", &TimeManager::setTimeStepSize);
}
}  // namespace vox