//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/wcsph/simulation_data_wcsph.h"
#include "vox.sph/wcsph/time_step_wcsph.h"

namespace py = pybind11;

void WCSPHModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Simulation Data WCSPH
    // ---------------------------------------
    py::class_<vox::SimulationDataWCSPH>(m_sub, "SimulationDataWCSPH")
            .def(py::init<>())
            .def("init", &vox::SimulationDataWCSPH::init)
            .def("cleanup", &vox::SimulationDataWCSPH::cleanup)
            .def("reset", &vox::SimulationDataWCSPH::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataWCSPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::SimulationDataWCSPH::emittedParticles)
            .def("getPressure", (Real(vox::SimulationDataWCSPH::*)(const unsigned int, const unsigned int) const) &
                                        vox::SimulationDataWCSPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataWCSPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &vox::SimulationDataWCSPH::setPressure)
            .def("getPressureAccel",
                 (const Vector3r& (vox::SimulationDataWCSPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataWCSPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataWCSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataWCSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &vox::SimulationDataWCSPH::setPressureAccel);

    // ---------------------------------------
    // Time Step WCSPH
    // ---------------------------------------
    py::class_<vox::TimeStepWCSPH, vox::TimeStep>(m_sub, "TimeStepWCSPH")
            .def_readwrite_static("STIFFNESS", &vox::TimeStepWCSPH::STIFFNESS)
            .def_readwrite_static("EXPONENT", &vox::TimeStepWCSPH::EXPONENT)
            .def(py::init<>());
}