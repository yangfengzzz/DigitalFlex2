//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "py.sph/common.h"
#include "vox.sph/iisph/simulation_data_iisph.h"
#include "vox.sph/iisph/time_step_iisph.h"

namespace py = pybind11;

void IISPHModule(const py::module &m_sub) {
    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<vox::SimulationDataIISPH>(m_sub, "SimulationDataIISPH")
            .def(py::init<>())
            .def("init", &vox::SimulationDataIISPH::init)
            .def("cleanup", &vox::SimulationDataIISPH::cleanup)
            .def("reset", &vox::SimulationDataIISPH::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataIISPH::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::SimulationDataIISPH::emittedParticles)

            .def("getAii", (Real(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int)
                                    const)(&vox::SimulationDataIISPH::getAii))
            // .def("getAii", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataIISPH::getAii)) // TODO: wont work by reference
            .def("setAii", &vox::SimulationDataIISPH::setAii)

            .def("getDii",
                 (const Vector3r &(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataIISPH::getDii)
            // .def("getDii", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getDii) // TODO: wont work by reference
            .def("setDii", &vox::SimulationDataIISPH::setDii)

            .def("getDij_pj",
                 (const Vector3r &(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataIISPH::getDij_pj)
            // .def("getDij_pj", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getDij_pj) // TODO: wont work by reference
            .def("setDij_pj", &vox::SimulationDataIISPH::setDij_pj)

            .def("getDensityAdv", (Real(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                                          vox::SimulationDataIISPH::getDensityAdv)
            // .def("getDensityAdv", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getDensityAdv) // TODO: wont work by reference
            .def("setDensityAdv", &vox::SimulationDataIISPH::setDensityAdv)

            .def("getPressure", (Real(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                                        vox::SimulationDataIISPH::getPressure)
            // .def("getPressure", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getPressure) // TODO: wont work by reference
            .def("setPressure", &vox::SimulationDataIISPH::setPressure)

            .def("getLastPressure", (Real(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                                            vox::SimulationDataIISPH::getLastPressure)
            // .def("getLastPressure", (Real& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getLastPressure) // TODO: wont work by reference
            .def("setLastPressure", &vox::SimulationDataIISPH::setLastPressure)

            .def("getPressureAccel",
                 (const Vector3r &(vox::SimulationDataIISPH::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataIISPH::getPressureAccel)
            // .def("getPressureAccel", (Vector3r& (SPH::SimulationDataIISPH::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataIISPH::getPressureAccel) // TODO: wont work by reference
            .def("setPressureAccel", &vox::SimulationDataIISPH::setPressureAccel);

    // ---------------------------------------
    // Class Simulation Data IISPH
    // ---------------------------------------
    py::class_<vox::TimeStepIISPH, vox::TimeStep>(m_sub, "TimeStepIISPH")
            .def("getSimulationData", &vox::TimeStepIISPH::getSimulationData)
            .def(py::init<>());
}
