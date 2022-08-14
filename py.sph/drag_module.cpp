//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/drag/drag_base.h"
#include "vox.sph/drag/drag_force_gissler2017.h"
#include "vox.sph/drag/drag_force_macklin2014.h"

namespace py = pybind11;

void DragModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Drag Base
    // ---------------------------------------
    py::class_<vox::DragBase, vox::NonPressureForceBase>(m_sub, "DragBase")
            .def_readwrite_static("DRAG_COEFFICIENT", &vox::DragBase::DRAG_COEFFICIENT);

    // ---------------------------------------
    // Class Drag Force Gissler
    // ---------------------------------------
    py::class_<vox::DragForce_Gissler2017, vox::DragBase>(m_sub, "DragForce_Gissler2017")
            .def(py::init<vox::FluidModel*>());

    // ---------------------------------------
    // Class Drag Force Macklin
    // ---------------------------------------
    py::class_<vox::DragForce_Macklin2014, vox::DragBase>(m_sub, "DragForce_Macklin2014")
            .def(py::init<vox::FluidModel*>());
}
