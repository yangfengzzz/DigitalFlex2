//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>

#include "py.sph/common.h"
#include "vox.sph/boundary_model.h"
#include "vox.sph/boundary_model_akinci2012.h"
#include "vox.sph/boundary_model_bender2019.h"
#include "vox.sph/boundary_model_koschier2017.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void BoundaryModelModule(const py::module& m_sub) {
    // ---------------------------------------
    // Boundary Model
    // ---------------------------------------
    py::class_<vox::BoundaryModel>(m_sub, "BoundaryModel")
            .def(py::init<>())
            .def("reset", &vox::BoundaryModel::reset)
            .def("performNeighborhoodSearchSort", &vox::BoundaryModel::performNeighborhoodSearchSort)
            .def("saveState", &vox::BoundaryModel::saveState)
            .def("loadState", &vox::BoundaryModel::loadState)
            .def("getRigidBodyObject", &vox::BoundaryModel::getRigidBodyObject,
                 py::return_value_policy::reference_internal)
            .def("addForce", overload_cast_<const Vector3r&, const Vector3r&>()(&vox::BoundaryModel::addForce))
#ifdef USE_AVX
            .def("addForce", overload_cast_<const Vector3f8&, const Vector3f8&, const unsigned int>()(
                                     &vox::BoundaryModel::addForce))
#endif
            .def("getPointVelocity", &vox::BoundaryModel::getPointVelocity)    // TODO: remove or fix
            .def("getForceAndTorque", &vox::BoundaryModel::getForceAndTorque)  // TODO: remove or fix
            .def("clearForceAndTorque", &vox::BoundaryModel::clearForceAndTorque);

    // ---------------------------------------
    // Boundary Model Akinci 2012
    // ---------------------------------------
    py::class_<vox::BoundaryModel_Akinci2012, vox::BoundaryModel>(m_sub, "BoundaryModelAkinci2012")
            .def(py::init<>())
            .def("numberOfParticles", &vox::BoundaryModel_Akinci2012::numberOfParticles)
            .def("getPointSetIndex", &vox::BoundaryModel_Akinci2012::getPointSetIndex)
            .def("computeBoundaryVolume", &vox::BoundaryModel_Akinci2012::computeBoundaryVolume)
            .def("resize", &vox::BoundaryModel_Akinci2012::resize)
            .def("reset", &vox::BoundaryModel_Akinci2012::reset)
            .def("performNeighborhoodSearchSort", &vox::BoundaryModel_Akinci2012::performNeighborhoodSearchSort)
            .def("saveState", &vox::BoundaryModel_Akinci2012::saveState)
            .def("loadState", &vox::BoundaryModel_Akinci2012::loadState)
            .def("initModel", &vox::BoundaryModel_Akinci2012::initModel)
            // .def("getPosition0", (Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned
            // int))(&vox::BoundaryModel_Akinci2012::getPosition0)) // TODO: wont work by reference
            .def("getPosition0", (const Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned int)
                                          const)(&vox::BoundaryModel_Akinci2012::getPosition0))
            .def("setPosition0", &vox::BoundaryModel_Akinci2012::setPosition0)
            // .def("getPosition", (Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned
            // int))(&vox::BoundaryModel_Akinci2012::getPosition)) TODO: wont work by reference
            .def("getPosition", (const Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned int)
                                         const)(&vox::BoundaryModel_Akinci2012::getPosition))
            .def("setPosition", &vox::BoundaryModel_Akinci2012::setPosition)
            // .def("getVelocity", (Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned
            // int))(&vox::BoundaryModel_Akinci2012::getVelocity)) // TODO: wont work by reference
            .def("getVelocity", (const Vector3r& (vox::BoundaryModel_Akinci2012::*)(const unsigned int)
                                         const)(&vox::BoundaryModel_Akinci2012::getVelocity))
            .def("setVelocity", &vox::BoundaryModel_Akinci2012::setVelocity)
            // .def("getVolume", (Real& (vox::BoundaryModel_Akinci2012::*)(const unsigned
            // int))(&vox::BoundaryModel_Akinci2012::getVolume)) // TODO: might work by reference, but not intended
            // behaviour. Use setter instead
            .def("getVolume", (const Real& (vox::BoundaryModel_Akinci2012::*)(const unsigned int)
                                       const)(&vox::BoundaryModel_Akinci2012::getVolume))
            .def("setVolume", &vox::BoundaryModel_Akinci2012::setVolume);

    // ---------------------------------------
    // Boundary Model Bender 2019
    // ---------------------------------------
    py::class_<vox::BoundaryModel_Bender2019, vox::BoundaryModel>(m_sub, "BoundaryModelBender2019")
            .def(py::init<>())
            .def("initModel", &vox::BoundaryModel_Bender2019::initModel)
            .def("reset", &vox::BoundaryModel_Bender2019::reset)
            .def("getMap", &vox::BoundaryModel_Bender2019::getMap, py::return_value_policy::reference_internal)
            .def("setMap", &vox::BoundaryModel_Bender2019::setMap)
            .def("getMaxDist", &vox::BoundaryModel_Bender2019::getMaxDist)
            .def("setMaxDist", &vox::BoundaryModel_Bender2019::setMaxDist)
            .def("getMaxVel", &vox::BoundaryModel_Bender2019::getMaxVel)
            .def("setMaxVel", &vox::BoundaryModel_Bender2019::setMaxVel)
            .def("getBoundaryVolume",
                 (const Real& (vox::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int)
                          const)(&vox::BoundaryModel_Bender2019::getBoundaryVolume))
            // .def("getBoundaryVolume", (Real& (vox::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned
            // int))(&vox::BoundaryModel_Bender2019::getBoundaryVolume)) // TODO: might work by reference, but not
            // intended behaviour. Use setter instead
            .def("setBoundaryVolume", &vox::BoundaryModel_Bender2019::setBoundaryVolume)
            .def("getBoundaryXj",
                 (const Vector3r& (vox::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned int)
                          const)(&vox::BoundaryModel_Bender2019::getBoundaryXj))
            // .def("getBoundaryXj", (Vector3r& (vox::BoundaryModel_Bender2019::*)(const unsigned int, const unsigned
            // int))(&vox::BoundaryModel_Bender2019::getBoundaryXj)) // TODO: wont work by reference
            .def("setBoundaryXj", &vox::BoundaryModel_Bender2019::setBoundaryXj);

    // ---------------------------------------
    // Boundary Model Koschier 2017
    // ---------------------------------------
    py::class_<vox::BoundaryModel_Koschier2017, vox::BoundaryModel>(m_sub, "BoundaryModelKoschier2017")
            .def(py::init<>())
            .def("initModel", &vox::BoundaryModel_Koschier2017::initModel)
            .def("reset", &vox::BoundaryModel_Koschier2017::reset)
            .def("getMap", &vox::BoundaryModel_Koschier2017::getMap, py::return_value_policy::reference_internal)
            .def("setMap", &vox::BoundaryModel_Koschier2017::setMap)
            .def("getMaxDist", &vox::BoundaryModel_Koschier2017::getMaxDist)
            .def("setMaxDist", &vox::BoundaryModel_Koschier2017::setMaxDist)
            .def("getMaxVel", &vox::BoundaryModel_Koschier2017::getMaxVel)
            .def("setMaxVel", &vox::BoundaryModel_Koschier2017::setMaxVel)
            .def("getBoundaryDensity",
                 (const Real& (vox::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)
                          const)(&vox::BoundaryModel_Koschier2017::getBoundaryDensity))
            // .def("getBoundaryDensity", (Real & (vox::BoundaryModel_Koschier2017::*)(const unsigned int, const
            // unsigned int))(&vox::BoundaryModel_Koschier2017::getBoundaryDensity)) // TODO: wont work by reference
            .def("setBoundaryDensity", &vox::BoundaryModel_Koschier2017::setBoundaryDensity)
            .def("getBoundaryDensityGradient",
                 (const Vector3r& (vox::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)
                          const)(&vox::BoundaryModel_Koschier2017::getBoundaryDensityGradient))
            // .def("getBoundaryDensityGradient", (Vector3r& (vox::BoundaryModel_Koschier2017::*)(const unsigned int,
            // const unsigned int))(&vox::BoundaryModel_Koschier2017::getBoundaryDensityGradient)) // TODO: wont work by
            // reference
            .def("setBoundaryDensityGradient", &vox::BoundaryModel_Koschier2017::setBoundaryDensityGradient)
            .def("getBoundaryXj",
                 (const Vector3r& (vox::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned int)
                          const)(&vox::BoundaryModel_Koschier2017::getBoundaryXj))
            // .def("getBoundaryXj", (Vector3r& (vox::BoundaryModel_Koschier2017::*)(const unsigned int, const unsigned
            // int))(&vox::BoundaryModel_Koschier2017::getBoundaryXj)) // TODO: wont work by reference
            .def("setBoundaryXj", &vox::BoundaryModel_Koschier2017::setBoundaryXj);
}
