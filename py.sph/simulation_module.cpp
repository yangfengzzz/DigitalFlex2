//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <string>
#include <vector>

#include "py.sph/common.h"
#include "vox.editor/boundary_simulator.h"
#include "vox.editor/exporter/exporter_base.h"
#include "vox.editor/gui/simulator_gui_base.h"
#include "vox.editor/static_boundary_simulator.h"
#include "vox.sph/scene_configuration.h"
#include "vox.sph/simulation.h"

namespace py = pybind11;
using namespace pybind11::literals;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void py_init_simulator(
        vox::SimulatorBase& obj,
        const std::string& sceneFile = "data/Scenes/DoubleDamBreak.json",  // TODO: change to empty default
        const std::string& programName = "pySPlisHSPlasH",
        bool useCache = true,
        const std::string& stateFile = "",
        const std::string& outputDir = "",
        bool initialPause = true,
        bool useGui = true,
        float stopAt = -1.0f,
        const std::string& param = "") {
    std::vector<std::string> argv;
    argv.push_back(programName);
    argv.emplace_back("--scene-file");
    argv.push_back(sceneFile);
    if (!useCache) argv.emplace_back("--no-cache");
    argv.emplace_back("--stopAt");
    argv.push_back(std::to_string(stopAt));
    if (strcmp(param.c_str(), "") != 0) {
        argv.emplace_back("--param");
        argv.push_back(param);
    }
    if (strcmp(stateFile.c_str(), "") != 0) {
        argv.emplace_back("--state-file");
        argv.push_back(stateFile);
    }
    if (strcmp(outputDir.c_str(), "") != 0) {
        argv.emplace_back("--output-dir");
        argv.push_back(outputDir);
    }
    if (!initialPause) argv.emplace_back("--no-initial-pause");
    if (!useGui) argv.emplace_back("--no-gui");
    obj.init(argv, "pySPlisHSPlasH");
}

void SimulationModule(py::module m_sub) {
    // ---------------------------------------
    // Enum Simulation Methods
    // ---------------------------------------
    py::enum_<vox::SimulationMethods>(m_sub, "SimulationMethods")
            .value("WCSPH", vox::SimulationMethods::WCSPH)
            .value("PCISPH", vox::SimulationMethods::PCISPH)
            .value("PBF", vox::SimulationMethods::PBF)
            .value("IISPH", vox::SimulationMethods::IISPH)
            .value("DFSPH", vox::SimulationMethods::DFSPH)
            .value("PF", vox::SimulationMethods::PF)
            .value("NumSimulationMethods", vox::SimulationMethods::NumSimulationMethods);

    // ---------------------------------------
    // Enum Boundary Handling Methods
    // ---------------------------------------
    py::enum_<vox::BoundaryHandlingMethods>(m_sub, "BoundaryHandlingMethods")
            .value("Akinci2012", vox::BoundaryHandlingMethods::Akinci2012)
            .value("Koschier2017", vox::BoundaryHandlingMethods::Koschier2017)
            .value("Bender2019", vox::BoundaryHandlingMethods::Bender2019)
            .value("NumSimulationMethods", vox::BoundaryHandlingMethods::NumSimulationMethods);

    // ---------------------------------------
    // Class Simulation
    // ---------------------------------------
    py::class_<vox::Simulation, vox::ParameterObject>(m_sub, "Simulation")
            .def_readwrite_static("SIM_2D", &vox::Simulation::SIM_2D)
            .def_readwrite_static("PARTICLE_RADIUS", &vox::Simulation::PARTICLE_RADIUS)
            .def_readwrite_static("GRAVITATION", &vox::Simulation::GRAVITATION)
            .def_readwrite_static("CFL_METHOD", &vox::Simulation::CFL_METHOD)
            .def_readwrite_static("CFL_FACTOR", &vox::Simulation::CFL_FACTOR)
            .def_readwrite_static("CFL_MIN_TIMESTEPSIZE", &vox::Simulation::CFL_MIN_TIMESTEPSIZE)
            .def_readwrite_static("CFL_MAX_TIMESTEPSIZE", &vox::Simulation::CFL_MAX_TIMESTEPSIZE)
            .def_readwrite_static("ENABLE_Z_SORT", &vox::Simulation::ENABLE_Z_SORT)

            .def_readwrite_static("KERNEL_METHOD", &vox::Simulation::KERNEL_METHOD)
            .def_readwrite_static("GRAD_KERNEL_METHOD", &vox::Simulation::GRAD_KERNEL_METHOD)
            .def_readwrite_static("ENUM_KERNEL_CUBIC", &vox::Simulation::ENUM_KERNEL_CUBIC)
            .def_readwrite_static("ENUM_KERNEL_WENDLANDQUINTICC2", &vox::Simulation::ENUM_KERNEL_WENDLANDQUINTICC2)
            .def_readwrite_static("ENUM_KERNEL_POLY6", &vox::Simulation::ENUM_KERNEL_POLY6)
            .def_readwrite_static("ENUM_KERNEL_SPIKY", &vox::Simulation::ENUM_KERNEL_SPIKY)
            .def_readwrite_static("ENUM_KERNEL_PRECOMPUTED_CUBIC", &vox::Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC)
            .def_readwrite_static("ENUM_KERNEL_CUBIC_2D", &vox::Simulation::ENUM_KERNEL_CUBIC_2D)
            .def_readwrite_static("ENUM_KERNEL_WENDLANDQUINTICC2_2D",
                                  &vox::Simulation::ENUM_KERNEL_WENDLANDQUINTICC2_2D)
            .def_readwrite_static("ENUM_GRADKERNEL_CUBIC", &vox::Simulation::ENUM_GRADKERNEL_CUBIC)
            .def_readwrite_static("ENUM_GRADKERNEL_WENDLANDQUINTICC2",
                                  &vox::Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2)
            .def_readwrite_static("ENUM_GRADKERNEL_POLY6", &vox::Simulation::ENUM_GRADKERNEL_POLY6)
            .def_readwrite_static("ENUM_GRADKERNEL_SPIKY", &vox::Simulation::ENUM_GRADKERNEL_SPIKY)
            .def_readwrite_static("ENUM_GRADKERNEL_PRECOMPUTED_CUBIC",
                                  &vox::Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC)
            .def_readwrite_static("ENUM_GRADKERNEL_CUBIC_2D", &vox::Simulation::ENUM_GRADKERNEL_CUBIC_2D)
            .def_readwrite_static("ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D",
                                  &vox::Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D)

            .def_readwrite_static("SIMULATION_METHOD", &vox::Simulation::SIMULATION_METHOD)

            .def_readwrite_static("ENUM_CFL_NONE", &vox::Simulation::ENUM_CFL_NONE)
            .def_readwrite_static("ENUM_CFL_STANDARD", &vox::Simulation::ENUM_CFL_STANDARD)
            .def_readwrite_static("ENUM_CFL_ITER", &vox::Simulation::ENUM_CFL_ITER)

            .def_readwrite_static("ENUM_SIMULATION_WCSPH", &vox::Simulation::ENUM_SIMULATION_WCSPH)
            .def_readwrite_static("ENUM_SIMULATION_PCISPH", &vox::Simulation::ENUM_SIMULATION_PCISPH)
            .def_readwrite_static("ENUM_SIMULATION_PBF", &vox::Simulation::ENUM_SIMULATION_PBF)
            .def_readwrite_static("ENUM_SIMULATION_IISPH", &vox::Simulation::ENUM_SIMULATION_IISPH)
            .def_readwrite_static("ENUM_SIMULATION_DFSPH", &vox::Simulation::ENUM_SIMULATION_DFSPH)
            .def_readwrite_static("ENUM_SIMULATION_PF", &vox::Simulation::ENUM_SIMULATION_PF)

            .def_readwrite_static("BOUNDARY_HANDLING_METHOD", &vox::Simulation::BOUNDARY_HANDLING_METHOD)
            .def_readwrite_static("ENUM_AKINCI2012", &vox::Simulation::ENUM_AKINCI2012)
            .def_readwrite_static("ENUM_KOSCHIER2017", &vox::Simulation::ENUM_KOSCHIER2017)
            .def_readwrite_static("ENUM_BENDER2019", &vox::Simulation::ENUM_BENDER2019)

            .def(py::init<>())
            .def("init", &vox::Simulation::init)
            .def("reset", &vox::Simulation::reset)
            .def_static("getCurrent", &vox::Simulation::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &vox::Simulation::setCurrent)
            .def_static("hasCurrent", &vox::Simulation::hasCurrent)

            .def("addFluidModel",
                 [](vox::Simulation& current, const std::string& id, std::vector<Vector3r> fluidParticles,
                    std::vector<Vector3r> fluidVelocities, std::vector<unsigned int> fluidObjectIds,
                    const unsigned int nMaxEmitterParticles) {
                     if (fluidParticles.size() != fluidVelocities.size())
                         throw std::runtime_error("Sizes of position and velocity array must be equal");
                     if (fluidParticles.size() != fluidObjectIds.size())
                         throw std::runtime_error("Sizes of position and object id array must be equal");
                     current.addFluidModel(id, fluidParticles.size(), fluidParticles.data(), fluidVelocities.data(),
                                           fluidObjectIds.data(), nMaxEmitterParticles);
                 })
            .def("getFluidModel", &vox::Simulation::getFluidModel, py::return_value_policy::reference_internal)
            .def("getFluidModelFromPointSet", &vox::Simulation::getFluidModelFromPointSet,
                 py::return_value_policy::reference_internal)
            .def("numberOfFluidModels", &vox::Simulation::numberOfFluidModels)

            .def("addBoundaryModel", &vox::Simulation::addBoundaryModel)
            .def("getBoundaryModel", &vox::Simulation::getBoundaryModel, py::return_value_policy::reference_internal)
            .def("getBoundaryModelFromPointSet", &vox::Simulation::getBoundaryModelFromPointSet,
                 py::return_value_policy::reference_internal)
            .def("numberOfBoundaryModels", &vox::Simulation::numberOfBoundaryModels)
            .def("updateBoundaryVolume", &vox::Simulation::updateBoundaryVolume)

            .def("getAnimationFieldSystem", &vox::Simulation::getAnimationFieldSystem,
                 py::return_value_policy::reference_internal)

            .def("getBoundaryHandlingMethod", &vox::Simulation::getBoundaryHandlingMethod)
            .def("setBoundaryHandlingMethod", &vox::Simulation::setBoundaryHandlingMethod)

            .def("getKernel", &vox::Simulation::getKernel)
            .def("setKernel", &vox::Simulation::setKernel)
            .def("getGradKernel", &vox::Simulation::getGradKernel)
            .def("setGradKernel", &vox::Simulation::setGradKernel)

            .def("W_zero", &vox::Simulation::W_zero)
            .def("W", &vox::Simulation::W)
            .def("gradW", &vox::Simulation::gradW)

            .def("getSimulationMethod", &vox::Simulation::getSimulationMethod)
            .def("setSimulationMethod", &vox::Simulation::setSimulationMethod)
            .def("setSimulationMethodChangedCallback", &vox::Simulation::setSimulationMethodChangedCallback)
            .def("getTimeStep", &vox::Simulation::getTimeStep,
                 py::return_value_policy::reference_internal)  // TODO: This returns abstract class pointer, figure out
                                                               // what to do with it
            .def("is2DSimulation", &vox::Simulation::is2DSimulation)
            .def("zSortEnabled", &vox::Simulation::zSortEnabled)
            .def("setParticleRadius", &vox::Simulation::setParticleRadius)
            .def("getParticleRadius", &vox::Simulation::getParticleRadius)
            .def("getSupportRadius", &vox::Simulation::getSupportRadius)
            .def("updateTimeStepSize", &vox::Simulation::updateTimeStepSize)
            .def("updateTimeStepSizeCFL", &vox::Simulation::updateTimeStepSizeCFL)
            .def("performNeighborhoodSearch", &vox::Simulation::performNeighborhoodSearch)
            .def("performNeighborhoodSearchSort", &vox::Simulation::performNeighborhoodSearchSort)
            .def("computeNonPressureForces", &vox::Simulation::computeNonPressureForces)
            .def("animateParticles", &vox::Simulation::animateParticles)
            .def("emitParticles", &vox::Simulation::emitParticles)
            .def("emittedParticles", &vox::Simulation::emittedParticles)
            .def("getNeighborhoodSearch", &vox::Simulation::getNeighborhoodSearch,
                 py::return_value_policy::reference_internal)
            .def("saveState", &vox::Simulation::saveState)
            .def("loadState", &vox::Simulation::loadState)
            .def("addViscosityMethod", &vox::Simulation::addViscosityMethod)
            .def("getViscosityMethods", &vox::Simulation::getViscosityMethods)
            .def("addVorticityMethod", &vox::Simulation::addVorticityMethod)
            .def("getVorticityMethods", &vox::Simulation::getVorticityMethods)
            .def("addDragMethod", &vox::Simulation::addDragMethod)
            .def("getDragMethods", &vox::Simulation::getDragMethods)
            .def("numberOfPointSets", &vox::Simulation::numberOfPointSets)
            .def("numberOfNeighbors", &vox::Simulation::numberOfNeighbors)
            .def("getNeighbor", &vox::Simulation::getNeighbor)
            .def("getNeighborList",
                 &vox::Simulation::getNeighborList);  // TODO: Might not work because of array pointer

    // ---------------------------------------
    // EXEC SUBMODULE
    // ---------------------------------------
    m_sub = m_sub.def_submodule("Exec");

    // ---------------------------------------
    // Struct Exporter
    // ---------------------------------------
    py::class_<vox::SimulatorBase::Exporter>(m_sub, "Exporter")
            .def_readwrite("key", &vox::SimulatorBase::Exporter::m_key)
            .def_readwrite("name", &vox::SimulatorBase::Exporter::m_name)
            .def_readwrite("decription", &vox::SimulatorBase::Exporter::m_description)
            .def_readwrite("id", &vox::SimulatorBase::Exporter::m_id)
            .def_readwrite("exporter", &vox::SimulatorBase::Exporter::m_exporter);

    // ---------------------------------------
    // Simulator Base class
    // ---------------------------------------
    py::class_<vox::SimulatorBase::SimulationMethod>(m_sub, "SimulationMethod")
            .def_readwrite("simulationMethod", &vox::SimulatorBase::SimulationMethod::simulationMethod)
            .def_readwrite("simulation", &vox::SimulatorBase::SimulationMethod::simulation)
            .def_readonly("model",
                          &vox::SimulatorBase::SimulationMethod::model);  // TODO: this is public property but defined
                                                                          // as readonly because of deleted assignment
                                                                          // operator. check in future

    py::class_<vox::SimulatorBase, vox::ParameterObject>(m_sub, "SimulatorBase")
            .def_readwrite_static("PAUSE", &vox::SimulatorBase::PAUSE)
            .def_readwrite_static("PAUSE_AT", &vox::SimulatorBase::PAUSE_AT)
            .def_readwrite_static("STOP_AT", &vox::SimulatorBase::STOP_AT)
            .def_readwrite_static("NUM_STEPS_PER_RENDER", &vox::SimulatorBase::NUM_STEPS_PER_RENDER)
            .def_readwrite_static("DATA_EXPORT_FPS", &vox::SimulatorBase::DATA_EXPORT_FPS)
            .def_readwrite_static("PARTICLE_EXPORT_ATTRIBUTES", &vox::SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES)
            .def_readwrite_static("STATE_EXPORT", &vox::SimulatorBase::STATE_EXPORT)
            .def_readwrite_static("STATE_EXPORT_FPS", &vox::SimulatorBase::STATE_EXPORT_FPS)
            .def_readwrite_static("RENDER_WALLS", &vox::SimulatorBase::RENDER_WALLS)

            .def_readwrite_static("ENUM_WALLS_NONE", &vox::SimulatorBase::ENUM_WALLS_NONE)
            .def_readwrite_static("ENUM_WALLS_PARTICLES_ALL", &vox::SimulatorBase::ENUM_WALLS_PARTICLES_ALL)
            .def_readwrite_static("ENUM_WALLS_PARTICLES_NO_WALLS", &vox::SimulatorBase::ENUM_WALLS_PARTICLES_NO_WALLS)
            .def_readwrite_static("ENUM_WALLS_GEOMETRY_ALL", &vox::SimulatorBase::ENUM_WALLS_GEOMETRY_ALL)
            .def_readwrite_static("ENUM_WALLS_GEOMETRY_NO_WALLS", &vox::SimulatorBase::ENUM_WALLS_GEOMETRY_NO_WALLS)

            .def(py::init<>())
            .def("run", &vox::SimulatorBase::run)
            // .def("init", [](vox::SimulatorBase& obj, std::vector<std::string> argv, std::string windowName){
            //     std::vector<const char *> cargv;
            //     cargv.reserve(argv.size());
            //     for (auto & elem : argv) {
            //         cargv.push_back(elem.c_str());
            //     }
            //     obj.init(argv.size(), const_cast<char**>(cargv.data()), windowName);
            // })
            .def("init", overload_cast_<std::vector<std::string>, const std::string&>()(&vox::SimulatorBase::init))
            .def("init", &py_init_simulator, "sceneFile"_a = "data/Scenes/DoubleDamBreak.json",
                 "programName"_a = "pySPlisHSPlasH", "useCache"_a = true, "stateFile"_a = "", "outputDir"_a = "",
                 "initialPause"_a = true, "useGui"_a = true, "stopAt"_a = -1.0, "param"_a = "")
            .def("initSimulation", &vox::SimulatorBase::initSimulation)
            .def("runSimulation", &vox::SimulatorBase::runSimulation)
            .def("cleanup", &vox::SimulatorBase::cleanup)

            .def("reset", &vox::SimulatorBase::reset)
            .def("timeStep", &vox::SimulatorBase::timeStep)
            .def("timeStepNoGUI", &vox::SimulatorBase::timeStepNoGUI)

            .def_static("particleInfo", &vox::SimulatorBase::particleInfo)
            .def("initDensityMap", &vox::SimulatorBase::initDensityMap)
            .def("initVolumeMap", &vox::SimulatorBase::initVolumeMap)

            .def("readParameters", &vox::SimulatorBase::readParameters)
            .def("step", &vox::SimulatorBase::step)

            .def("saveState", &vox::SimulatorBase::saveState)
            // .def("loadStateDialog", &vox::SimulatorBase::loadStateDialog)
            .def("loadState", &vox::SimulatorBase::loadState)
            .def("writeFluidParticlesState", &vox::SimulatorBase::writeFluidParticlesState)
            .def("readFluidParticlesState", &vox::SimulatorBase::readFluidParticlesState)
            .def("writeBoundaryState", &vox::SimulatorBase::writeBoundaryState)
            .def("readBoundaryState", &vox::SimulatorBase::readBoundaryState)
            .def("writeParameterState", &vox::SimulatorBase::writeParameterState)
            .def("readParameterState", &vox::SimulatorBase::readParameterState)
            .def("writeParameterObjectState", &vox::SimulatorBase::writeParameterObjectState)
            .def("readParameterObjectState", &vox::SimulatorBase::readParameterObjectState)

            .def("updateBoundaryParticles", &vox::SimulatorBase::updateBoundaryParticles)
            .def("updateDMVelocity", &vox::SimulatorBase::updateDMVelocity)
            .def("updateVMVelocity", &vox::SimulatorBase::updateVMVelocity)

            .def("getScalarField", &vox::SimulatorBase::getScalarField)
            .def("updateScalarField", &vox::SimulatorBase::updateScalarField)
            .def("determineMinMaxOfScalarField", &vox::SimulatorBase::determineMinMaxOfScalarField)

            .def_static("loadObj", &vox::SimulatorBase::loadObj)

            .def("getSceneLoader", &vox::SimulatorBase::getSceneLoader, py::return_value_policy::reference_internal)

            .def("getExePath", &vox::SimulatorBase::getExePath)

            .def("getUseParticleCaching", &vox::SimulatorBase::getUseParticleCaching)
            .def("setUseParticleCaching", &vox::SimulatorBase::setUseParticleCaching)
            .def("getUseGUI", &vox::SimulatorBase::getUseGUI)
            .def("setUseGUI", &vox::SimulatorBase::setUseGUI)

            .def("getColorField", &vox::SimulatorBase::getColorField)
            .def("setColorField", &vox::SimulatorBase::setColorField)

            .def("getColorMapType", &vox::SimulatorBase::getColorMapType)
            .def("setColorMapType", &vox::SimulatorBase::setColorMapType)
            .def("getRenderMaxValue", &vox::SimulatorBase::getRenderMaxValue)
            .def("setRenderMaxValue", &vox::SimulatorBase::setRenderMaxValue)
            .def("getRenderMinValue", &vox::SimulatorBase::getRenderMinValue)
            .def("setRenderMinValue", &vox::SimulatorBase::setRenderMinValue)
            .def("getOutputPath", &vox::SimulatorBase::getOutputPath)

            .def("getStateFile", &vox::SimulatorBase::getStateFile)
            .def("setStateFile", &vox::SimulatorBase::setStateFile)

            .def("getBoundarySimulator", &vox::SimulatorBase::getBoundarySimulator,
                 py::return_value_policy::reference_internal)
            .def("setBoundarySimulator", &vox::SimulatorBase::setBoundarySimulator)
            .def("getGui", &vox::SimulatorBase::getGui, py::return_value_policy::reference_internal)
            .def("setGui", &vox::SimulatorBase::setGui)
            .def("isStaticScene", &vox::SimulatorBase::isStaticScene)

            .def("addParticleExporter", &vox::SimulatorBase::addParticleExporter)
            .def("getParticleExporters", &vox::SimulatorBase::getParticleExporters)
            .def("addRigidBodyExporter", &vox::SimulatorBase::addRigidBodyExporter)
            .def("getRigidBodyExporters", &vox::SimulatorBase::getRigidBodyExporters)

            .def("activateExporter", &vox::SimulatorBase::activateExporter)

            .def("setTimeStepCB", &vox::SimulatorBase::setTimeStepCB)
            .def("setResetCB", &vox::SimulatorBase::setResetCB);

    // ---------------------------------------
    // SceneConfiguration
    // ---------------------------------------
    py::class_<vox::SceneConfiguration>(m_sub, "SceneConfiguration")
            .def_static("getCurrent", &vox::SceneConfiguration::getCurrent, py::return_value_policy::reference)
            .def_static("setCurrent", &vox::SceneConfiguration::setCurrent)
            .def_static("hasCurrent", &vox::SceneConfiguration::hasCurrent)
            .def("getSceneFile", &vox::SceneConfiguration::getSceneFile)
            .def("getScene", &vox::SceneConfiguration::getScene, py::return_value_policy::reference_internal);

    // ---------------------------------------
    // BoundarySimulator
    // ---------------------------------------
    py::class_<vox::BoundarySimulator>(m_sub, "BoundarySimulator")
            .def(py::init<>())
            .def("init", &vox::BoundarySimulator::init)
            .def("timeStep", &vox::BoundarySimulator::timeStep)
            .def("initBoundaryData", &vox::BoundarySimulator::initBoundaryData)
            .def("reset", &vox::BoundarySimulator::reset)
            .def("updateBoundaryForces", &vox::BoundarySimulator::updateBoundaryForces);

    // ---------------------------------------
    // StaticBoundarySimulator
    // ---------------------------------------
    py::class_<vox::StaticBoundarySimulator, vox::BoundarySimulator>(m_sub, "StaticBoundarySimulator")
            .def(py::init<vox::SimulatorBase*>());
}
