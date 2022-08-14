//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/common.h"
#include "vox.sph/gui/opengl/simulator_opengl.h"
#include "vox.sph/position_based_dynamics_wrapper/pbd_boundary_simulator.h"
#include "vox.sph/simulator_base.h"
#ifdef USE_IMGUI
#include "vox.sph/gui/imgui/simulator_gui_imgui.h.h"
#include "vox.sph/position_based_dynamics_wrapper/pbd_simulator_gui_imgui.h.h"
#else
#include "vox.sph/gui/tweakbar/simulator_gui_tweakbar.h"
#include "vox.sph/position_based_dynamics_wrapper/pbd_simulator_gui_tweakbar.h"
#endif

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW
#endif
#endif

using namespace vox;
using namespace std;

SimulatorBase *base = nullptr;
Simulator_GUI_Base *gui = nullptr;

// main
int main(int argc, char **argv) {
    REPORT_MEMORY_LEAKS;
    base = new SimulatorBase();
    base->init(argc, argv, "SPlisHSPlasH");

    if (base->getUseGUI()) {
#ifdef USE_IMGUI
        if (base->isStaticScene())
            gui = new Simulator_GUI_imgui(base);
        else
            gui = new PBD_Simulator_GUI_imgui(base,
                                              ((PBDBoundarySimulator *)base->getBoundarySimulator())->getPBDWrapper());
#else
        if (base->isStaticScene())
            gui = new Simulator_GUI_TweakBar(base);
        else
            gui = new PBD_Simulator_GUI_TweakBar(
                    base, ((PBDBoundarySimulator *)base->getBoundarySimulator())->getPBDWrapper());
#endif
        base->setGui(gui);
    }
    base->run();

    delete base;
    delete gui;

    return 0;
}
