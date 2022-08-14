//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/projective_fluids/simulation_data_pf.h"
#include "vox.sph/projective_fluids/time_step_pf.h"

namespace py = pybind11;

void PFModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Simulation Data PF
    // ---------------------------------------
    py::class_<vox::SimulationDataPF>(m_sub, "SimulationDataPF")
            .def(py::init<>())
            .def("init", &vox::SimulationDataPF::init)
            .def("cleanup", &vox::SimulationDataPF::cleanup)
            .def("reset", &vox::SimulationDataPF::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataPF::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::SimulationDataPF::emittedParticles)

            .def("getOldPosition", (Vector3r(vox::SimulationDataPF::*)(const unsigned int, const unsigned int)
                                            const)(&vox::SimulationDataPF::getOldPosition))
            // .def("getOldPosition", (Vector3r& (vox::SimulationDataPF::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataPF::getOldPosition)) // TODO: wont work by reference
            .def("setOldPosition", &vox::SimulationDataPF::setOldPosition)

            .def("getNumFluidNeighbors",
                 (unsigned int (vox::SimulationDataPF::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPF::getNumFluidNeighbors)
            // .def("getNumFluidNeighbors", (unsigned int& (vox::SimulationDataPF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPF::getNumFluidNeighbors) // TODO: wont work by reference
            .def("setNumFluidNeighbors", &vox::SimulationDataPF::setNumFluidNeighbors)

            .def("getS", (const Vector3r& (vox::SimulationDataPF::*)(const unsigned int, const unsigned int) const) &
                                 vox::SimulationDataPF::getS)
            // .def("getS", (Vector3r& (vox::SimulationDataPF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPF::getS) // TODO: wont work by reference
            .def("setS", &vox::SimulationDataPF::setS)

            .def("getDiag", (const Vector3r& (vox::SimulationDataPF::*)(const unsigned int, const unsigned int) const) &
                                    vox::SimulationDataPF::getDiag)
            // .def("getDiag", (Vector3r& (vox::SimulationDataPF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPF::getDiag) // TODO: wont work by reference
            .def("setDiag", &vox::SimulationDataPF::setDiag)

            .def("getParticleOffset", &vox::SimulationDataPF::getParticleOffset);

    // ---------------------------------------
    // Class Time Step PF
    // ---------------------------------------
    py::class_<vox::TimeStepPF, vox::TimeStep>(m_sub, "TimeStepPF")
            .def_readwrite_static("STIFFNESS", &vox::TimeStepPF::STIFFNESS)
            .def(py::init<>())
            .def_static("matrixVecProd", &vox::TimeStepPF::matrixVecProd);
}