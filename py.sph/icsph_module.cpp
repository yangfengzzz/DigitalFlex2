//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "py.sph/common.h"
#include "vox.sph/icsph/simulation_data_icsph.h"
#include "vox.sph/icsph/time_step_icsph.h"

namespace py = pybind11;

void ICSPHModule(const py::module &m_sub) {
    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<vox::SimulationDataICSPH>(m_sub, "SimulationDataICSPH")
            .def(py::init<>())
            .def("init", &vox::SimulationDataICSPH::init)
            .def("cleanup", &vox::SimulationDataICSPH::cleanup)
            .def("reset", &vox::SimulationDataICSPH::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataICSPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::SimulationDataICSPH::emittedParticles)

            .def("getAii", (Real(vox::SimulationDataICSPH::*)(const unsigned int, const unsigned int)
                                    const)(&vox::SimulationDataICSPH::getAii))
            // .def("getAii", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataICSPH::getAii)) // TODO: wont work by reference
            .def("setAii", &vox::SimulationDataICSPH::setAii)

            .def("getDensityAdv", (Real(vox::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const) &
                                          vox::SimulationDataICSPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataICSPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &vox::SimulationDataICSPH::setDensityAdv)

            .def("getPressure", (Real(vox::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const) &
                                        vox::SimulationDataICSPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataICSPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &vox::SimulationDataICSPH::setPressure)

            .def("getPressureAccel",
                 (const Vector3r &(vox::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataICSPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataICSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &vox::SimulationDataICSPH::setPressureAccel)

            .def("getPressureGradient",
                 (const Vector3r &(vox::SimulationDataICSPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataICSPH::getPressureGradient)
            // .def("getPressureGradient", (Vector3r& (SPH::SimulationDataICSPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataICSPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureGradient", &vox::SimulationDataICSPH::setPressureGradient);

    // ---------------------------------------
    // Class Simulation Data ICSPH
    // ---------------------------------------
    py::class_<vox::TimeStepICSPH, vox::TimeStep>(m_sub, "TimeStepICSPH")
            .def_readwrite_static("LAMBDA", &vox::TimeStepICSPH::LAMBDA)
            .def_readwrite_static("PRESSURE_CLAMPING", &vox::TimeStepICSPH::PRESSURE_CLAMPING)
            .def("getSimulationData", &vox::TimeStepICSPH::getSimulationData)
            .def(py::init<>());
}
