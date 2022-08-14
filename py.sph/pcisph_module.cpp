//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/pcisph/simulation_data_pcisph.h"
#include "vox.sph/pcisph/time_step_pcisph.h"

namespace py = pybind11;

void PCISPHModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Simulation Data PCISPH
    // ---------------------------------------
    py::class_<vox::SimulationDataPCISPH>(m_sub, "SimulationDataPCISPH")
            .def(py::init<>())
            .def("init", &vox::SimulationDataPCISPH::init)
            .def("cleanup", &vox::SimulationDataPCISPH::cleanup)
            .def("reset", &vox::SimulationDataPCISPH::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataPCISPH::performNeighborhoodSearchSort)
            .def("getPCISPH_ScalingFactor", &vox::SimulationDataPCISPH::getPCISPH_ScalingFactor)
            .def("emittedParticles", &vox::SimulationDataPCISPH::emittedParticles)

            .def("getPredictedPosition",
                 (const Vector3r& (vox::SimulationDataPCISPH::*)(const unsigned int, const unsigned int)
                          const)(&vox::SimulationDataPCISPH::getPredictedPosition))
            // .def("getPredictedPosition", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPCISPH::getPredictedPosition) // TODO: wont work by reference
            .def("setPredictedPosition", &vox::SimulationDataPCISPH::setPredictedPosition)

            .def("getPredictedVelocity",
                 (const Vector3r& (vox::SimulationDataPCISPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPCISPH::getPredictedVelocity)
            // .def("getPredictedVelocity", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPCISPH::getPredictedVelocity) // TODO: wont work by reference
            .def("setPredictedVelocity", &vox::SimulationDataPCISPH::setPredictedVelocity)

            .def("getDensityAdv", (Real(vox::SimulationDataPCISPH::*)(const unsigned int, const unsigned int) const) &
                                          vox::SimulationDataPCISPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPCISPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &vox::SimulationDataPCISPH::setDensityAdv)

            .def("getPressure", (Real(vox::SimulationDataPCISPH::*)(const unsigned int, const unsigned int) const) &
                                        vox::SimulationDataPCISPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPCISPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &vox::SimulationDataPCISPH::setPressure)

            .def("getPressureAccel",
                 (const Vector3r& (vox::SimulationDataPCISPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPCISPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataPCISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPCISPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &vox::SimulationDataPCISPH::setPressureAccel);

    // ---------------------------------------
    // Class Time Step PCISPH
    // ---------------------------------------
    py::class_<vox::TimeStepPCISPH, vox::TimeStep>(m_sub, "TimeStepPCISPH").def(py::init<>());
}