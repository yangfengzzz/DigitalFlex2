//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include "py.sph/common.h"
#include "vox.sph/elasticity/elasticity_base.h"
#include "vox.sph/elasticity/elasticity_becker2009.h"
#include "vox.sph/elasticity/elasticity_kugelstadt2021.h"
#include "vox.sph/elasticity/elasticity_peer2018.h"

namespace py = pybind11;

void ElasticityModule(const py::module& m_sub) {
    // ---------------------------------------
    // Class Elasticity Base
    // ---------------------------------------
    py::class_<vox::ElasticityBase, vox::NonPressureForceBase>(m_sub, "ElasticityBase")
            .def_readwrite_static("YOUNGS_MODULUS", &vox::ElasticityBase::YOUNGS_MODULUS)
            .def_readwrite_static("POISSON_RATIO", &vox::ElasticityBase::POISSON_RATIO)
            .def_readwrite_static("POISSON_RATIO", &vox::ElasticityBase::FIXED_BOX_MIN)
            .def_readwrite_static("POISSON_RATIO", &vox::ElasticityBase::FIXED_BOX_MAX);

    // ---------------------------------------
    // Class Elasticity Becker 2009
    // ---------------------------------------
    py::class_<vox::Elasticity_Becker2009, vox::ElasticityBase>(m_sub, "Elasticity_Becker2009")
            .def_readwrite_static("ALPHA", &vox::Elasticity_Becker2009::ALPHA)
            .def(py::init<vox::FluidModel*>());

    // ---------------------------------------
    // Class Elasticity Peer 2018
    // ---------------------------------------
    py::class_<vox::Elasticity_Peer2018, vox::ElasticityBase>(m_sub, "Elasticity_Peer2018")
            .def_readwrite_static("ITERATIONS", &vox::Elasticity_Peer2018::ITERATIONS)
            .def_readwrite_static("MAX_ITERATIONS", &vox::Elasticity_Peer2018::MAX_ITERATIONS)
            .def_readwrite_static("MAX_ERROR", &vox::Elasticity_Peer2018::MAX_ERROR)
            .def_readwrite_static("ALPHA", &vox::Elasticity_Peer2018::ALPHA)

            .def_static("matrixVecProd", &vox::Elasticity_Peer2018::matrixVecProd)
            .def(py::init<vox::FluidModel*>());

    py::class_<vox::Elasticity_Kugelstadt2021, vox::ElasticityBase>(m_sub, "Elasticity_Kugelstadt2021")
            .def_readwrite_static("ITERATIONS_V", &vox::Elasticity_Kugelstadt2021::ITERATIONS_V)
            .def_readwrite_static("MAX_ITERATIONS_V", &vox::Elasticity_Kugelstadt2021::MAX_ITERATIONS_V)
            .def_readwrite_static("MAX_ERROR_V", &vox::Elasticity_Kugelstadt2021::MAX_ERROR_V)
            .def_readwrite_static("ALPHA", &vox::Elasticity_Kugelstadt2021::ALPHA)
            .def_readwrite_static("MAX_NEIGHBORS", &vox::Elasticity_Kugelstadt2021::MAX_NEIGHBORS)

            .def_static("matrixVecProd", &vox::Elasticity_Kugelstadt2021::matrixVecProd)
            .def("computeRotations", &vox::Elasticity_Kugelstadt2021::computeRotations)
            .def(py::init<vox::FluidModel*>());
}
