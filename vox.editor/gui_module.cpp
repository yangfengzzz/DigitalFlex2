//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "py.sph/common.h"
#include "vox.editor/gui/simulator_gui_base.h"
#ifdef USE_IMGUI
#include "vox.editor/gui/imgui/simulator_gui_imgui.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_simulator_gui_imgui.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_wrapper.h"
#else
#include "vox.editor/gui/tweakbar/simulator_gui_tweakbar.h"
#include "vox.editor/position_based_dynamics_wrapper/pbd_simulator_gui_tweakbar.h"
#endif
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>

#include "vox.editor/gui/opengl/simulator_opengl.h"

namespace py = pybind11;

void GUIModule(py::module m) {
    // ---------------------------------------
    // GUI Submodule
    // ---------------------------------------
    auto m_sub = m.def_submodule("GUI");

    // ---------------------------------------
    // Simulator GUI Base Class
    // ---------------------------------------
    py::class_<vox::Simulator_GUI_Base>(m_sub, "Simulator_GUI_Base")
            .def(py::init<vox::SimulatorBase*>())
            .def("init",
                 [](vox::Simulator_GUI_Base& obj, std::vector<std::string> argv, std::string windowName) {
                     std::vector<const char*> cargv;
                     cargv.reserve(argv.size());
                     for (auto& elem : argv) cargv.push_back(elem.c_str());
                     obj.init(static_cast<int>(argv.size()), const_cast<char**>(cargv.data()), windowName.c_str());
                 })
            .def("initSimulationParameterGUI", &vox::Simulator_GUI_Base::initSimulationParameterGUI)
            .def("initParameterGUI", &vox::Simulator_GUI_Base::initParameterGUI)
            .def("render", &vox::Simulator_GUI_Base::render)
            .def("reset", &vox::Simulator_GUI_Base::reset)
            .def("update", &vox::Simulator_GUI_Base::update)
            .def("cleanup", &vox::Simulator_GUI_Base::cleanup)
            .def("run", &vox::Simulator_GUI_Base::run)
            .def("stop", &vox::Simulator_GUI_Base::stop)
            .def("addKeyFunc", &vox::Simulator_GUI_Base::addKeyFunc)
            .def("getSimulatorBase", &vox::Simulator_GUI_Base::getSimulatorBase,
                 py::return_value_policy::reference_internal);

#ifdef USE_IMGUI
    // ---------------------------------------
    // imgui GUI Class
    // ---------------------------------------
    py::class_<vox::Simulator_GUI_imgui, vox::Simulator_GUI_Base>(m_sub, "Simulator_GUI_imgui")
            .def(py::init<vox::SimulatorBase*>());

    // ---------------------------------------
    // pbd imgui GUI Class
    // ---------------------------------------
    py::class_<vox::PBD_Simulator_GUI_imgui, vox::Simulator_GUI_imgui>(m_sub, "PBD_Simulator_GUI_imgui")
            .def(py::init<vox::SimulatorBase*, PBDWrapper*>());

#else
    // ---------------------------------------
    // Tweakbar GUI Class
    // ---------------------------------------
    py::class_<vox::Simulator_GUI_TweakBar, vox::Simulator_GUI_Base>(m_sub, "Simulator_GUI_TweakBar")
            .def(py::init<vox::SimulatorBase*>());
    // .def("getTweakBar", &vox::Simulator_GUI_TweakBar::getTweakBar); //TODO: Wrap tweak bar for this

    // ---------------------------------------
    // pbd tweakbar GUI Class
    // ---------------------------------------
    py::class_<vox::PBD_Simulator_GUI_TweakBar, vox::Simulator_GUI_TweakBar>(m_sub, "PBD_Simulator_GUI_TweakBar")
            .def(py::init<vox::SimulatorBase*, vox::PBDWrapper*>());
#endif

    // ---------------------------------------
    // OpenGL GUI Class // TODO implement missing functions
    // ---------------------------------------
    py::class_<vox::Simulator_OpenGL>(m_sub, "Simulator_OpenGL")
            .def(py::init<>())
            .def("initShaders", &vox::Simulator_OpenGL::initShaders)
            .def("renderFluid", &vox::Simulator_OpenGL::renderFluid)
            .def("renderSelectedParticles", &vox::Simulator_OpenGL::renderSelectedParticles)
            .def("renderBoundary", &vox::Simulator_OpenGL::renderBoundary)
            .def("renderBoundaryParticles", &vox::Simulator_OpenGL::renderBoundaryParticles);
}
