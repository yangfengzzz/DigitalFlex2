//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

// Put this here for now
#include "py.pbd/common.h"
#include "vox.base/logging.h"
#include "vox.base/timing.h"

namespace py = pybind11;

void CollisionDetectionModule(py::module);
void ConstraintsModule(py::module);
void ParameterObjectModule(py::module);
void ParticleDataModule(py::module);
void RigidBodyModule(py::module);
void SimulationModelModule(py::module);
void SimulationModule(py::module);
void TimeStepModule(py::module);
void TimeModule(py::module);
void UtilitiesModule(py::module);

PYBIND11_MODULE(pypbd, m) {
    CollisionDetectionModule(m);
    ParticleDataModule(m);
    RigidBodyModule(m);
    ParameterObjectModule(m);
    SimulationModelModule(m);
    SimulationModule(m);
    ConstraintsModule(m);
    TimeStepModule(m);
    TimeModule(m);
    UtilitiesModule(m);
}