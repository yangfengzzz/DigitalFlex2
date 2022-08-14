//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include <utility>

#include "py.pbd/common.h"
#include "vox.pbd/cubic_sdf_collision_detection.h"
#include "vox.pbd/simulation.h"

namespace py = pybind11;

void SimulationModule(py::module m_sub) {
    // ---------------------------------------
    // Class Simulation
    // ---------------------------------------
    py::class_<vox::Simulation, vox::ParameterObject>(std::move(m_sub), "Simulation")
            .def(py::init<>())
            .def_static("getCurrent", &vox::Simulation::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &vox::Simulation::setCurrent)
            .def_static("hasCurrent", &vox::Simulation::hasCurrent)

            .def("init", &vox::Simulation::init)
            .def("reset", &vox::Simulation::reset)
            .def("getModel", &vox::Simulation::getModel, py::return_value_policy::reference_internal)
            .def("setModel", &vox::Simulation::setModel)
            .def("getTimeStep", &vox::Simulation::getTimeStep, py::return_value_policy::reference_internal)
            .def("setTimeStep", &vox::Simulation::setTimeStep)
            .def("initDefault", [](vox::Simulation& sim) {
                sim.setModel(new vox::SimulationModel());
                sim.getModel()->init();
                auto* cd = new vox::CubicSDFCollisionDetection();
                sim.getTimeStep()->setCollisionDetection(*sim.getModel(), cd);
            });
}