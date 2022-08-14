//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/pbf/simulation_data_pbf.h"
#include "vox.sph/pbf/time_integration.h"
#include "vox.sph/pbf/time_step_pbf.h"

namespace py = pybind11;

void PBFModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Simulation Data PBF
    // ---------------------------------------
    py::class_<vox::SimulationDataPBF>(m_sub, "SimulationDataPBF")
            .def(py::init<>())
            .def("init", &vox::SimulationDataPBF::init)
            .def("cleanup", &vox::SimulationDataPBF::cleanup)
            .def("reset", &vox::SimulationDataPBF::reset)
            .def("performNeighborhoodSearchSort", &vox::SimulationDataPBF::performNeighborhoodSearchSort)
            .def("emittedParticles", &vox::SimulationDataPBF::emittedParticles)

            .def("getLambda", (const Real& (vox::SimulationDataPBF::*)(const unsigned int, const unsigned int)
                                       const)(&vox::SimulationDataPBF::getLambda))
            // .def("getLambda", (Real & (vox::SimulationDataPBF::*)(const unsigned int, const unsigned
            // int))(&vox::SimulationDataPBF::getLambda)) // TODO: wont work by reference
            .def("setLambda", &vox::SimulationDataPBF::setLambda)

            .def("getDeltaX",
                 (const Vector3r& (vox::SimulationDataPBF::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPBF::getDeltaX)
            // .def("getDeltaX", (Vector3r & (vox::SimulationDataPBF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPBF::getDeltaX) // TODO: wont work by reference
            .def("setDeltaX", &vox::SimulationDataPBF::setDeltaX)

            .def("getLastPosition",
                 (const Vector3r& (vox::SimulationDataPBF::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPBF::getLastPosition)
            // .def("getLastPosition", (Vector3r & (vox::SimulationDataPBF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPBF::getLastPosition) // TODO: wont work by reference
            .def("setLastPosition", &vox::SimulationDataPBF::setLastPosition)

            .def("getOldPosition",
                 (const Vector3r& (vox::SimulationDataPBF::*)(const unsigned int, const unsigned int) const) &
                         vox::SimulationDataPBF::getOldPosition)
            // .def("getOldPosition", (Vector3r & (vox::SimulationDataPBF::*)(const unsigned int, const unsigned
            // int))&vox::SimulationDataPBF::getOldPosition) // TODO: wont work by reference
            .def("setOldPosition", &vox::SimulationDataPBF::setOldPosition);

    // ---------------------------------------
    // Class Time Integration
    // ---------------------------------------
    py::class_<vox::TimeIntegration>(m_sub, "TimeIntegration")
            .def_static("semiImplicitEuler", &vox::TimeIntegration::semiImplicitEuler)
            .def_static("velocityUpdateFirstOrder", &vox::TimeIntegration::velocityUpdateFirstOrder)
            .def_static("velocityUpdateSecondOrder", &vox::TimeIntegration::velocityUpdateSecondOrder);

    // ---------------------------------------
    // Class Time Step PBF
    // ---------------------------------------
    py::class_<vox::TimeStepPBF, vox::TimeStep>(m_sub, "TimeStepPBF")
            .def_readwrite_static("VELOCITY_UPDATE_METHOD", &vox::TimeStepPBF::VELOCITY_UPDATE_METHOD)
            .def_readwrite_static("ENUM_PBF_FIRST_ORDER", &vox::TimeStepPBF::ENUM_PBF_FIRST_ORDER)
            .def_readwrite_static("ENUM_PBF_SECOND_ORDER", &vox::TimeStepPBF::ENUM_PBF_SECOND_ORDER)
            .def(py::init<>());
}