//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/dfsph/simulation_data_dfsph.h"
#include "vox.sph/dfsph/time_step_dfsph.h"

namespace py = pybind11;

void DFSPHModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Simulation Data DFSPH
    // ---------------------------------------
    py::class_<vox::SimulationDataDFSPH>(m_sub, "SimulationDataDFSPH")
            .def(py::init<>())
            .def("init", &vox::SimulationDataDFSPH::init)
            .def("cleanup", &vox::SimulationDataDFSPH::cleanup)
            .def("reset", &vox::SimulationDataDFSPH::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataDFSPH::performNeighborhoodSearchSort)

            .def("getFactor", (Real(vox::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)
                                       const)(&vox::SimulationDataDFSPH::getFactor))
            // .def("getFactor", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataDFSPH::getFactor)) // TODO: wont work by reference
            .def("setFactor", &vox::SimulationDataDFSPH::setFactor)

            .def("getKappa", (Real(vox::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)
                                      const)(&vox::SimulationDataDFSPH::getKappa))
            // .def("getKappa", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataDFSPH::getKappa)) // TODO: wont work by reference
            .def("setKappa", &vox::SimulationDataDFSPH::setKappa)

            .def("getKappaV", (Real(vox::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)
                                       const)(&vox::SimulationDataDFSPH::getKappaV))
            // .def("getKappaV", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataDFSPH::getKappaV)) // TODO: wont work by reference
            .def("setKappaV", &vox::SimulationDataDFSPH::setKappaV)

            .def("getDensityAdv", (Real(vox::SimulationDataDFSPH::*)(const unsigned int, const unsigned int)
                                           const)(&vox::SimulationDataDFSPH::getDensityAdv))
            // .def("getDensityAdv", (Real& (SPH::SimulationDataDFSPH::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataDFSPH::getDensityAdv)) // TODO: wont work by reference
            .def("setDensityAdv", &vox::SimulationDataDFSPH::setDensityAdv);

    // ---------------------------------------
    // Class Time Step DFSPH
    // ---------------------------------------
    py::class_<vox::TimeStepDFSPH, vox::TimeStep>(m_sub, "TimeStepDFSPH")
            .def_readwrite_static("SOLVER_ITERATIONS_V", &vox::TimeStepDFSPH::SOLVER_ITERATIONS_V)
            .def_readwrite_static("MAX_ITERATIONS_V", &vox::TimeStepDFSPH::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &vox::TimeStepDFSPH::MAX_ERROR_V)
            .def_readwrite_static("USE_DIVERGENCE_SOLVER", &vox::TimeStepDFSPH::USE_DIVERGENCE_SOLVER)

            .def(py::init<>());
}