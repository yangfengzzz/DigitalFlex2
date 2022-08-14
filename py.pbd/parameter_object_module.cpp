//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>

#include <utility>

#include "vox.base/parameter_object.h"
#include "vox.pbd/common.h"

namespace py = pybind11;

void ParameterObjectModule(py::module m_sub) {
    // auto m_sub = m.def_submodule("Common");
    py::class_<vox::ParameterObject>(std::move(m_sub), "ParameterObject")
            // .def(py::init<>()) TODO: no constructor for now because this object does not need to be constructable
            .def("getValueBool", &vox::ParameterObject::getValue<bool>)
            .def("getValueInt", &vox::ParameterObject::getValue<int>)
            .def("getValueUInt", &vox::ParameterObject::getValue<unsigned int>)
            .def("getValueFloat", &vox::ParameterObject::getValue<Real>)
            .def("getValueString", &vox::ParameterObject::getValue<std::string>)

            .def("setValueBool", &vox::ParameterObject::setValue<bool>)
            .def("setValueInt", &vox::ParameterObject::setValue<int>)
            .def("setValueUInt", &vox::ParameterObject::setValue<unsigned int>)
            .def("setValueFloat", &vox::ParameterObject::setValue<Real>)
            .def("setValueString", &vox::ParameterObject::setValue<std::string>);
}
