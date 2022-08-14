//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "py.sph/bind_pointer_vector.h"
#include "vox.sph/animation_field.h"
#include "vox.sph/animation_field_system.h"

namespace py = pybind11;

void AnimationFieldModule(const py::module &m) {
    py::class_<vox::AnimationField>(m, "AnimationField")
            .def(py::init<>([](const std::string &particleFieldName, const Vector3r &pos, const Matrix3r &rotation,
                               const Vector3r &scale, const std::string expression[3], const unsigned int type = 0) {
                return vox::AnimationField(particleFieldName, pos, rotation, scale, expression, type);
            }))
            .def("setStartTime", &vox::AnimationField::setStartTime)
            .def("setEndTime", &vox::AnimationField::setEndTime)
            .def("step", &vox::AnimationField::step);

    py::bind_pointer_vector<std::vector<vox::AnimationField *>>(m, "AnimationFieldVector");

    py::class_<vox::AnimationFieldSystem>(m, "AnimationFieldSystem")
            .def(py::init<>())
            .def("addAnimationField", &vox::AnimationFieldSystem::addAnimationField)
            .def("numAnimationFields", &vox::AnimationFieldSystem::numAnimationFields)
            .def("getAnimationFields", &vox::AnimationFieldSystem::getAnimationFields,
                 py::return_value_policy::reference_internal)
            .def("step", &vox::AnimationFieldSystem::step)
            .def("reset", &vox::AnimationFieldSystem::reset);
}