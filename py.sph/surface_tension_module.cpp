//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/non_pressure_force_base.h"
#include "vox.sph/surface_tension/surface_tension_akinci2013.h"
#include "vox.sph/surface_tension/surface_tension_becker2007.h"
#include "vox.sph/surface_tension/surface_tension_he2014.h"

namespace py = pybind11;

void SurfaceTensionModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Surface Tension Base
    // ---------------------------------------
    py::class_<vox::SurfaceTensionBase, vox::NonPressureForceBase>(m_sub, "SurfaceTensionBase")
            .def_readwrite_static("SURFACE_TENSION", &vox::SurfaceTensionBase::SURFACE_TENSION)
            .def_readwrite_static("SURFACE_TENSION_BOUNDARY", &vox::SurfaceTensionBase::SURFACE_TENSION_BOUNDARY);

    // ---------------------------------------
    // Class Surface Akinci 2013
    // ---------------------------------------
    py::class_<vox::SurfaceTension_Akinci2013, vox::SurfaceTensionBase>(m_sub, "SurfaceTension_Akinci2013")
            .def(py::init<vox::FluidModel*>())
            .def("computeNormals", &vox::SurfaceTension_Akinci2013::computeNormals)
            // .def("getNormal", (Vector3r& (vox::SurfaceTension_Akinci2013::*)(const unsigned
            // int))&vox::SurfaceTension_Akinci2013::getNormal) // TODO: wont work by reference
            .def("getNormal", (const Vector3r& (vox::SurfaceTension_Akinci2013::*)(const unsigned int) const) &
                                      vox::SurfaceTension_Akinci2013::getNormal)
            .def("setNormal", &vox::SurfaceTension_Akinci2013::setNormal);

    // ---------------------------------------
    // Class Surface Becker 2007
    // ---------------------------------------
    py::class_<vox::SurfaceTension_Becker2007, vox::SurfaceTensionBase>(m_sub, "SurfaceTension_Becker2007")
            .def(py::init<vox::FluidModel*>());

    // ---------------------------------------
    // Class Surface He 2014
    // ---------------------------------------
    py::class_<vox::SurfaceTension_He2014, vox::SurfaceTensionBase>(m_sub, "SurfaceTension_He2014")
            .def(py::init<vox::FluidModel*>())
            .def("getColor",
                 (Real(vox::SurfaceTension_He2014::*)(const unsigned int) const)(&vox::SurfaceTension_He2014::getColor))
            // .def("getColor", (Real& (vox::SurfaceTension_He2014::*)(const unsigned
            // int))(&vox::SurfaceTension_He2014::getColor)) // TODO: wont work by reference
            .def("setColor", &vox::SurfaceTension_He2014::setColor)
            .def("getGradC2", (Real(vox::SurfaceTension_He2014::*)(const unsigned int)
                                       const)(&vox::SurfaceTension_He2014::getGradC2))
            // .def("getGradC2", (Real& (vox::SurfaceTension_He2014::*)(const unsigned
            // int))(&vox::SurfaceTension_He2014::getGradC2)) // TODO: wont work by reference
            .def("setGradC2", &vox::SurfaceTension_He2014::setGradC2);
}