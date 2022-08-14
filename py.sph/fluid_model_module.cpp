//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include "py.sph/common.h"
#include "vox.sph/drag/drag_base.h"
#include "vox.sph/elasticity/elasticity_base.h"
#include "vox.sph/emitter_system.h"
#include "vox.sph/surface_tension/surface_tension_base.h"
#include "vox.sph/viscosity/viscosity_base.h"
#include "vox.sph/vorticity/vorticity_base.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

// TODO: remove reference getters

std::function<void*(const unsigned int)> makeVoidPointerFct(
        py::array_t<Real, py::array::c_style | py::array::forcecast> arr) {
    return [arr](const unsigned int i) mutable -> void* { return arr.mutable_unchecked().mutable_data(i); };
}

void FluidModelModule(py::module m_sub) {
    // ---------------------------------------
    // Enum Field Type
    // ---------------------------------------
    py::enum_<vox::FieldType>(m_sub, "FieldType")
            .value("Scalar", vox::FieldType::Scalar)
            .value("Vector3", vox::FieldType::Vector3)
            .value("Vector6", vox::FieldType::Vector6)
            .value("Matrix3", vox::FieldType::Matrix3)
            .value("Matrix6", vox::FieldType::Matrix6)
            .value("UInt", vox::FieldType::UInt);

    // ---------------------------------------
    // Struct Field Description
    // ---------------------------------------
    py::class_<vox::FieldDescription>(m_sub, "FieldDescription")
            .def(py::init<>([](const std::string& n, const vox::FieldType& t,
                               const std::function<void*(const unsigned int)>& fct,
                               const bool s = false) { return vox::FieldDescription(n, t, fct, s); }))
            .def_readwrite("name", &vox::FieldDescription::name)
            .def_readwrite("type", &vox::FieldDescription::type)
            .def_readwrite("getFct", &vox::FieldDescription::getFct)
            .def_readwrite("storeData", &vox::FieldDescription::storeData);

    py::bind_vector<std::vector<vox::FieldDescription>>(m_sub, "FieldDescriptionVector");

    // ---------------------------------------
    // Enum class Particle State
    // ---------------------------------------
    py::enum_<vox::ParticleState>(m_sub, "ParticleState")
            .value("Active", vox::ParticleState::Active)
            .value("AnimatedByEmitter", vox::ParticleState::AnimatedByEmitter);

    // ---------------------------------------
    // Class Fluid Model
    // ---------------------------------------
    py::class_<vox::FluidModel, vox::ParameterObject>(m_sub, "FluidModel")
            .def_readwrite_static("NUM_PARTICLES", &vox::FluidModel::NUM_PARTICLES)
            .def_readwrite_static("NUM_REUSED_PARTICLES", &vox::FluidModel::NUM_REUSED_PARTICLES)
            .def_readwrite_static("DENSITY0", &vox::FluidModel::DENSITY0)

            .def_readwrite_static("DRAG_METHOD", &vox::FluidModel::DRAG_METHOD)
            .def_readwrite_static("SURFACE_TENSION_METHOD", &vox::FluidModel::SURFACE_TENSION_METHOD)
            .def_readwrite_static("VISCOSITY_METHOD", &vox::FluidModel::VISCOSITY_METHOD)
            .def_readwrite_static("VORTICITY_METHOD", &vox::FluidModel::VORTICITY_METHOD)
            .def_readwrite_static("ELASTICITY_METHOD", &vox::FluidModel::ELASTICITY_METHOD)

            .def(py::init<>())
            .def("init", &vox::FluidModel::init)
            .def("getId", &vox::FluidModel::getId)
            .def("getDensity0", &vox::FluidModel::getDensity0)
            .def("setDensity0", &vox::FluidModel::setDensity0)
            .def("getPointSetIndex", &vox::FluidModel::getPointSetIndex)
            .def("addField", &vox::FluidModel::addField)
            .def("getFields", &vox::FluidModel::getFields,
                 py::return_value_policy::reference_internal)  // TODO: Bind return vector?
            .def("getFieldBuffer",
                 [](vox::FluidModel& obj, const unsigned int i) -> py::memoryview {
                     auto& field = obj.getField(i);
                     void* base_ptr = field.getFct(0);
                     int num_particles = obj.numParticles();
                     switch (field.type) {
                         case vox::FieldType::Scalar:
                             return py::memoryview::from_buffer((Real*)base_ptr, {num_particles}, {sizeof(Real)});
                         case vox::FieldType::Vector3:
                             return py::memoryview::from_buffer((Real*)base_ptr, {num_particles, 3},
                                                                {sizeof(Real) * 3, sizeof(Real)});
                         case vox::FieldType::UInt:
                             return py::memoryview::from_buffer((unsigned int*)base_ptr, {num_particles},
                                                                {sizeof(unsigned int)});
                         default:
                             break;
                     }
                     return py::memoryview(py::buffer_info());
                 })
            .def("getFieldBuffer",
                 [](vox::FluidModel& obj, const std::string& name) -> py::memoryview {
                     auto& field = obj.getField(name);
                     void* base_ptr = field.getFct(0);
                     int num_particles = obj.numParticles();
                     switch (field.type) {
                         case vox::FieldType::Scalar:
                             return py::memoryview::from_buffer((Real*)base_ptr, {num_particles}, {sizeof(Real)});
                         case vox::FieldType::Vector3:
                             return py::memoryview::from_buffer((Real*)base_ptr, {num_particles, 3},
                                                                {sizeof(Real) * 3, sizeof(Real)});
                         case vox::FieldType::UInt:
                             return py::memoryview::from_buffer((unsigned int*)base_ptr, {num_particles},
                                                                {sizeof(unsigned int)});
                         default:
                             break;
                     }
                     return py::memoryview(py::buffer_info());
                 })
            .def("getField", overload_cast_<const unsigned int>()(&vox::FluidModel::getField))
            .def("getField", overload_cast_<const unsigned int>()(&vox::FluidModel::getField))
            .def("getField", overload_cast_<const std::string&>()(&vox::FluidModel::getField))
            .def("numberOfFields", &vox::FluidModel::numberOfFields)
            .def("removeFieldByName", &vox::FluidModel::removeFieldByName)
            .def("setNumActiveParticles", &vox::FluidModel::setNumActiveParticles)
            .def("numberOfParticles", &vox::FluidModel::numberOfParticles)
            .def("getEmitterSystem", &vox::FluidModel::getEmitterSystem, py::return_value_policy::reference_internal)
            .def("reset", &vox::FluidModel::reset)
            .def("performNeighborhoodSearchSort", &vox::FluidModel::performNeighborhoodSearchSort)
            .def("initModel", &vox::FluidModel::initModel)
            .def("numParticles", &vox::FluidModel::numParticles)
            .def("numActiveParticles", &vox::FluidModel::numActiveParticles)
            .def("getNumActiveParticles0", &vox::FluidModel::getNumActiveParticles0)
            .def("setNumActiveParticles0", &vox::FluidModel::setNumActiveParticles0)
            .def("emittedParticles", &vox::FluidModel::emittedParticles)

            .def("getSurfaceTensionMethod", &vox::FluidModel::getSurfaceTensionMethod)
            .def("setSurfaceTensionMethod",
                 overload_cast_<const unsigned int>()(&vox::FluidModel::setSurfaceTensionMethod))
            .def("setSurfaceTensionMethod",
                 overload_cast_<const std::string&>()(&vox::FluidModel::setSurfaceTensionMethod))
            .def("getViscosityMethod", &vox::FluidModel::getViscosityMethod)
            .def("setViscosityMethod", overload_cast_<const unsigned int>()(&vox::FluidModel::setViscosityMethod))
            .def("setViscosityMethod", overload_cast_<const std::string&>()(&vox::FluidModel::setViscosityMethod))
            .def("getVorticityMethod", &vox::FluidModel::getVorticityMethod)
            .def("setVorticityMethod", overload_cast_<const unsigned int>()(&vox::FluidModel::setVorticityMethod))
            .def("setVorticityMethod", overload_cast_<const std::string&>()(&vox::FluidModel::setVorticityMethod))
            .def("getDragMethod", &vox::FluidModel::getDragMethod)
            .def("setDragMethod", overload_cast_<const unsigned int>()(&vox::FluidModel::setDragMethod))
            .def("setDragMethod", overload_cast_<const std::string&>()(&vox::FluidModel::setDragMethod))
            .def("getElasticityMethod", &vox::FluidModel::getElasticityMethod)
            .def("setElasticityMethod", overload_cast_<const unsigned int>()(&vox::FluidModel::setElasticityMethod))
            .def("setElasticityMethod", overload_cast_<const std::string&>()(&vox::FluidModel::setElasticityMethod))

            .def("getSurfaceTensionBase", &vox::FluidModel::getSurfaceTensionBase,
                 py::return_value_policy::reference_internal)
            .def("getViscosityBase", &vox::FluidModel::getViscosityBase, py::return_value_policy::reference_internal)
            .def("getVorticityBase", &vox::FluidModel::getVorticityBase, py::return_value_policy::reference_internal)
            .def("getDragBase", &vox::FluidModel::getDragBase, py::return_value_policy::reference_internal)
            .def("getElasticityBase", &vox::FluidModel::getElasticityBase, py::return_value_policy::reference_internal)

            .def("setDragMethodChangedCallback", &vox::FluidModel::setDragMethodChangedCallback)
            .def("setSurfaceMethodChangedCallback", &vox::FluidModel::setSurfaceMethodChangedCallback)
            .def("setViscosityMethodChangedCallback", &vox::FluidModel::setViscosityMethodChangedCallback)
            .def("setVorticityMethodChangedCallback", &vox::FluidModel::setVorticityMethodChangedCallback)
            .def("setElasticityMethodChangedCallback", &vox::FluidModel::setElasticityMethodChangedCallback)

            .def("computeSurfaceTension", &vox::FluidModel::computeSurfaceTension)
            .def("computeViscosity", &vox::FluidModel::computeViscosity)
            .def("computeVorticity", &vox::FluidModel::computeVorticity)
            .def("computeDragForce", &vox::FluidModel::computeDragForce)
            .def("computeElasticity", &vox::FluidModel::computeElasticity)

            .def("saveState", &vox::FluidModel::saveState)
            .def("loadState", &vox::FluidModel::loadState)

            // .def("getPosition0", (Vector3r& (vox::FluidModel::*)(const unsigned
            // int))(&vox::FluidModel::getPosition0)) // TODO: wont work by reference
            .def("getPosition0",
                 (const Vector3r& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getPosition0))
            .def("setPosition0", &vox::FluidModel::setPosition0)

            // .def("getPosition", (Vector3r& (vox::FluidModel::*)(const unsigned int))(&vox::FluidModel::getPosition))
            // // TODO: wont work by reference
            .def("getPosition",
                 (const Vector3r& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getPosition))
            .def("setPosition", &vox::FluidModel::setPosition)

            // .def("getVelocity", (Vector3r& (vox::FluidModel::*)(const unsigned int))(&vox::FluidModel::getVelocity))
            // // TODO: wont work by reference
            .def("getVelocity",
                 (const Vector3r& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getVelocity))
            .def("setVelocity", &vox::FluidModel::setVelocity)

            // .def("getVelocity0", (Vector3r& (vox::FluidModel::*)(const unsigned
            // int))(&vox::FluidModel::getVelocity0)) // TODO: wont work by reference
            .def("getVelocity0",
                 (const Vector3r& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getVelocity0))
            .def("setVelocity0", &vox::FluidModel::setVelocity0)

            // .def("getAcceleration", (Vector3r& (vox::FluidModel::*)(const unsigned
            // int))(&vox::FluidModel::getAcceleration)) // TODO: wont work by reference
            .def("getAcceleration",
                 (const Vector3r& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getAcceleration))
            .def("setAcceleration", &vox::FluidModel::setAcceleration)

            // .def("getMass", (Real& (vox::FluidModel::*)(const unsigned int))(&vox::FluidModel::getMass)) // TODO:
            // wont work by reference
            .def("getMass", (Real(vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getMass))
            .def("setMass", &vox::FluidModel::setMass)

            // .def("getDensity", (Real& (vox::FluidModel::*)(const unsigned int))(&vox::FluidModel::getDensity)) //
            // TODO: wont work by reference
            .def("getDensity",
                 (const Real& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getDensity))
            .def("setMass", &vox::FluidModel::setMass)

            // .def("getParticleId", (unsigned int& (vox::FluidModel::*)(const unsigned
            // int))(&vox::FluidModel::getParticleId)) // TODO: wont work by reference
            .def("getParticleId",
                 (const unsigned int& (vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getParticleId))

            // .def("getParticleState", (vox::ParticleState& (vox::FluidModel::*)(const unsigned
            // int))(&vox::FluidModel::getParticleState)) // TODO: wont work by reference
            .def("getParticleState", (const vox::ParticleState& (vox::FluidModel::*)(const unsigned int)
                                              const)(&vox::FluidModel::getParticleState))
            .def("setParticleState", &vox::FluidModel::setParticleState)

            // .def("getVolume", (Real& (vox::FluidModel::*)(const unsigned int))(&vox::FluidModel::getVolume)) // TODO:
            // wont work by reference
            .def("getVolume", (Real(vox::FluidModel::*)(const unsigned int) const)(&vox::FluidModel::getVolume));

    m_sub.def("makeVoidPointerFct", &makeVoidPointerFct, py::arg().noconvert());
}
