//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "common.h"
#include "vox.sph/viscosity/viscosity_bender2017.h"
#include "vox.sph/viscosity/viscosity_peer2015.h"
#include "vox.sph/viscosity/viscosity_peer2016.h"
#include "vox.sph/viscosity/viscosity_standard.h"
#include "vox.sph/viscosity/viscosity_takahashi2015.h"
#include "vox.sph/viscosity/viscosity_weiler2018.h"
#include "vox.sph/viscosity/viscosity_xsph.h"

namespace py = pybind11;

void ViscosityModule(const py::module& m_sub) {
    // ---------------------------------------
    // Viscosity Base
    // ---------------------------------------
    py::class_<vox::ViscosityBase, vox::NonPressureForceBase>(m_sub, "ViscosityBase")
            .def_readwrite_static("VISCOSITY_COEFFICIENT", &vox::ViscosityBase::VISCOSITY_COEFFICIENT);

    // ---------------------------------------
    // Viscosity Bender 2017
    // ---------------------------------------
    py::class_<vox::Viscosity_Bender2017, vox::ViscosityBase>(m_sub, "Viscosity_Bender2017")
            .def_readwrite_static("ITERATIONS", &vox::Viscosity_Bender2017::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::Viscosity_Bender2017::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::Viscosity_Bender2017::MAX_ERROR)

            .def(py::init<vox::FluidModel*>())
            .def("computeTargetStrainRate", &vox::Viscosity_Bender2017::computeTargetStrainRate)
            .def("computeViscosityFactor", &vox::Viscosity_Bender2017::computeViscosityFactor)
            .def("viscoGradientMultTransposeRightOpt", &vox::Viscosity_Bender2017::viscoGradientMultTransposeRightOpt)
            .def("getTargetStrainRate", (const Vector6r& (vox::Viscosity_Bender2017::*)(const unsigned int)
                                                 const)(&vox::Viscosity_Bender2017::getTargetStrainRate))
            // .def("getTargetStrainRate", (Vector6r& (vox::Viscosity_Bender2017::*)(const unsigned
            // int))&vox::Viscosity_Bender2017::getTargetStrainRate) // TODO: wont work by reference
            .def("setTargetStrainRate", &vox::Viscosity_Bender2017::setTargetStrainRate)
            .def("getViscosityFactor", (const Matrix6r& (vox::Viscosity_Bender2017::*)(const unsigned int) const) &
                                               vox::Viscosity_Bender2017::getViscosityFactor)
            // .def("getViscosityFactor", (Matrix6r& (vox::Viscosity_Bender2017::*)(const unsigned
            // int))&vox::Viscosity_Bender2017::getViscosityFactor) // TODO: wont work by reference
            .def("setViscosityFactor", &vox::Viscosity_Bender2017::setViscosityFactor)
            .def("getViscosityLambda", (const Vector6r& (vox::Viscosity_Bender2017::*)(const unsigned int) const) &
                                               vox::Viscosity_Bender2017::getViscosityLambda)
            // .def("getViscosityLambda", (Vector6r& (vox::Viscosity_Bender2017::*)(const unsigned
            // int))&vox::Viscosity_Bender2017::getViscosityLambda) // TODO: wont work by reference
            .def("setViscosityLambda", &vox::Viscosity_Bender2017::setViscosityLambda);

    // ---------------------------------------
    // Viscosity Peer 2015
    // ---------------------------------------
    py::class_<vox::Viscosity_Peer2015, vox::ViscosityBase>(m_sub, "Viscosity_Peer2015")
            .def_readwrite_static("ITERATIONS", &vox::Viscosity_Peer2015::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::Viscosity_Peer2015::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::Viscosity_Peer2015::MAX_ERROR)
            .def(py::init<vox::FluidModel*>())
            .def_static("matrixVecProd", &vox::Viscosity_Peer2015::matrixVecProd)
            .def_static("diagonalMatrixElement", &vox::Viscosity_Peer2015::diagonalMatrixElement)
            .def("getTargetNablaV", (const Matrix3r& (vox::Viscosity_Peer2015::*)(const unsigned int)
                                             const)(&vox::Viscosity_Peer2015::getTargetNablaV))
            // .def("getTargetNablaV", (Matrix3r& (vox::Viscosity_Peer2015::*)(const unsigned
            // int))&vox::Viscosity_Peer2015::getTargetNablaV) // TODO: wont work by reference
            .def("setTargetNablaV", &vox::Viscosity_Peer2015::setTargetNablaV);

    // ---------------------------------------
    // Viscosity Peer2016
    // ---------------------------------------
    py::class_<vox::Viscosity_Peer2016, vox::ViscosityBase>(m_sub, "Viscosity_Peer2016")
            .def_readwrite_static("ITERATIONS_V", &vox::Viscosity_Peer2016::ITERATIONS_V)
            .def_readwrite_static("ITERATIONS_OMEGA", &vox::Viscosity_Peer2016::ITERATIONS_OMEGA)
            .def_readwrite_static("MAX_ITERATIONS_V", &vox::Viscosity_Peer2016::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &vox::Viscosity_Peer2016::MAX_ERROR_V)
            .def_readwrite_static("MAX_ITERATIONS_OMEGA", &vox::Viscosity_Peer2016::MAX_ITERATIONS_OMEGA)
            .def_readwrite_static("MAX_ERROR_OMEGA", &vox::Viscosity_Peer2016::MAX_ERROR_OMEGA)

            .def(py::init<vox::FluidModel*>())
            .def_static("matrixVecProdV", &vox::Viscosity_Peer2016::matrixVecProdV)
            .def_static("diagonalMatrixElementV", &vox::Viscosity_Peer2016::diagonalMatrixElementV)
            .def_static("matrixVecProdOmega", &vox::Viscosity_Peer2016::matrixVecProdOmega)
            .def_static("diagonalMatrixElementOmega", &vox::Viscosity_Peer2016::diagonalMatrixElementOmega)
            .def("getTargetNablaV", (const Matrix3r& (vox::Viscosity_Peer2016::*)(const unsigned int)
                                             const)(&vox::Viscosity_Peer2016::getTargetNablaV))
            // .def("getTargetNablaV", (Matrix3r& (vox::Viscosity_Peer2016::*)(const unsigned
            // int))&vox::Viscosity_Peer2016::getTargetNablaV) // TODO: wont work by reference
            .def("setTargetNablaV", &vox::Viscosity_Peer2016::setTargetNablaV)
            .def("getOmega", (const Vector3r& (vox::Viscosity_Peer2016::*)(const unsigned int) const) &
                                     vox::Viscosity_Peer2016::getOmega)
            // .def("getOmega", ( Vector3r& (vox::Viscosity_Peer2016::*)(const unsigned
            // int))&vox::Viscosity_Peer2016::getOmega) // TODO: wont work by reference
            .def("setOmega", &vox::Viscosity_Peer2016::setOmega);

    // ---------------------------------------
    // Viscosity Standard
    // ---------------------------------------
    py::class_<vox::Viscosity_Standard, vox::ViscosityBase>(m_sub, "Viscosity_Standard")
            .def_readwrite_static("VISCOSITY_COEFFICIENT_BOUNDARY", &vox::Viscosity_Standard::VISCOSITY_COEFFICIENT)
            .def(py::init<vox::FluidModel*>());

    // ---------------------------------------
    // Viscosity
    // ---------------------------------------
    py::class_<vox::Viscosity_Takahashi2015, vox::ViscosityBase>(m_sub, "Viscosity_Takahashi2015")
            .def_readwrite_static("ITERATIONS", &vox::Viscosity_Takahashi2015::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::Viscosity_Takahashi2015::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::Viscosity_Takahashi2015::MAX_ERROR)

            .def(py::init<vox::FluidModel*>())
            .def_static("matrixVecProd", &vox::Viscosity_Takahashi2015::matrixVecProd)
            .def("getViscousStress", (const Matrix3r& (vox::Viscosity_Takahashi2015::*)(const unsigned int) const) &
                                             vox::Viscosity_Takahashi2015::getViscousStress)
            // .def("getViscousStress", (Matrix3r& (vox::Viscosity_Takahashi2015::*)(const unsigned
            // int))&vox::Viscosity_Takahashi2015::getViscousStress) // TODO: wont work by reference
            .def("setViscousStress", &vox::Viscosity_Takahashi2015::setViscousStress)
            .def("getAccel", (const Vector3r& (vox::Viscosity_Takahashi2015::*)(const unsigned int) const) &
                                     vox::Viscosity_Takahashi2015::getAccel)
            // .def("getAccel", (Vector3r& (vox::Viscosity_Takahashi2015::*)(const unsigned
            // int))&vox::Viscosity_Takahashi2015::getAccel) // TODO: wont work by reference
            .def("setAccel", &vox::Viscosity_Takahashi2015::setAccel);

    // ---------------------------------------
    // Viscosity Weiler 2018
    // ---------------------------------------
    py::class_<vox::Viscosity_Weiler2018, vox::ViscosityBase>(m_sub, "Viscosity_Weiler2018")
            .def_readwrite_static("ITERATIONS", &vox::Viscosity_Weiler2018::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::Viscosity_Weiler2018::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::Viscosity_Weiler2018::MAX_ERROR)
            .def_readwrite_static("VISCOSITY_COEFFICIENT_BOUNDARY",
                                  &vox::Viscosity_Weiler2018::VISCOSITY_COEFFICIENT_BOUNDARY)

            .def(py::init<vox::FluidModel*>())
            .def_static("matrixVecProd", &vox::Viscosity_Weiler2018::matrixVecProd)
            .def("getVDiff", (const Vector3r& (vox::Viscosity_Weiler2018::*)(const unsigned int) const) &
                                     vox::Viscosity_Weiler2018::getVDiff)
            .def("setVDiff", &vox::Viscosity_Weiler2018::setVDiff);

    // ---------------------------------------
    // Viscosity XSPH
    // ---------------------------------------
    py::class_<vox::Viscosity_XSPH, vox::ViscosityBase>(m_sub, "Viscosity_XSPH")
            .def_readwrite_static("VISCOSITY_COEFFICIENT_BOUNDARY",
                                  &vox::Viscosity_XSPH::VISCOSITY_COEFFICIENT_BOUNDARY)

            .def(py::init<vox::FluidModel*>());
}