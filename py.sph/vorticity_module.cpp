//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "common.h"
#include "vox.base/common.h"
#include "vox.sph/non_pressure_force_base.h"
#include "vox.sph/viscosity/viscosity_bender2017.h"
#include "vox.sph/vorticity/micropolar_model_bender2017.h"
#include "vox.sph/vorticity/vorticity_confinement.h"

namespace py = pybind11;

void VorticityModule(const py::module& m_sub) {
    // ---------------------------------------
    // Vorticity Base
    // ---------------------------------------
    py::class_<vox::VorticityBase, vox::NonPressureForceBase>(m_sub, "VorticityBase")
            .def_readwrite_static("VORTICITY_COEFFICIENT", &vox::VorticityBase::VORTICITY_COEFFICIENT);

    // ---------------------------------------
    // Vorticity Bender 2017
    // ---------------------------------------
    py::class_<vox::MicropolarModel_Bender2017, vox::VorticityBase>(m_sub, "MicropolarModel_Bender2017")
            .def_readwrite_static("VISCOSITY_OMEGA", &vox::MicropolarModel_Bender2017::VISCOSITY_OMEGA)
            .def_readwrite_static("INERTIA_INVERSE", &vox::MicropolarModel_Bender2017::INERTIA_INVERSE)

            .def(py::init<vox::FluidModel*>())
            .def("getAngularAcceleration",
                 (const Vector3r& (vox::MicropolarModel_Bender2017::*)(const unsigned int) const) &
                         vox::MicropolarModel_Bender2017::getAngularAcceleration)
            // .def("getAngularAcceleration", (Vector3r& (vox::MicropolarModel_Bender2017::*)(const unsigned
            // int))&vox::MicropolarModel_Bender2017::getAngularAcceleration) // TODO: wont work by reference
            .def("setAngularAcceleration", &vox::MicropolarModel_Bender2017::setAngularAcceleration)
            .def("getAngularVelocity",
                 (const Vector3r& (vox::MicropolarModel_Bender2017::*)(const unsigned int) const) &
                         vox::MicropolarModel_Bender2017::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (vox::MicropolarModel_Bender2017::*)(const unsigned
            // int))&vox::MicropolarModel_Bender2017::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &vox::MicropolarModel_Bender2017::setAngularVelocity);

    // ---------------------------------------
    // Vorticity Bender 2017
    // ---------------------------------------
    py::class_<vox::VorticityConfinement, vox::VorticityBase>(m_sub, "VorticityConfinement")
            .def(py::init<vox::FluidModel*>())
            .def("getAngularVelocity", (const Vector3r& (vox::VorticityConfinement::*)(const unsigned int) const) &
                                               vox::VorticityConfinement::getAngularVelocity)
            // .def("getAngularVelocity", (Vector3r& (vox::VorticityConfinement::*)(const unsigned
            // int))&vox::VorticityConfinement::getAngularVelocity) // TODO: wont work by reference
            .def("setAngularVelocity", &vox::VorticityConfinement::setAngularVelocity);
}