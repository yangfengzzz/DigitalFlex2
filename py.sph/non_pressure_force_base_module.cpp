//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include <utility>

#include "py.sph/common.h"
#include "vox.sph/non_pressure_force_base.h"

namespace py = pybind11;

void NonPressureForceBaseModule(py::module m_sub) {
    py::class_<vox::NonPressureForceBase, vox::ParameterObject>(std::move(m_sub), "NonPressureForceBase")
            .def("step", &vox::NonPressureForceBase::step)
            .def("reset", &vox::NonPressureForceBase::reset)
            .def("performNeighborhoodSearchSort", &vox::NonPressureForceBase::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::NonPressureForceBase::emittedParticles)
            .def("saveState", &vox::NonPressureForceBase::saveState)
            .def("loadState", &vox::NonPressureForceBase::loadState)
            .def("getModel", &vox::NonPressureForceBase::getModel, py::return_value_policy::reference_internal)
            .def("init", &vox::NonPressureForceBase::init);
}
