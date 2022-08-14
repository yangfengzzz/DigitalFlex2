//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/parameter_object.h"
#include "vox.editor/simulator_base.h"
#include "vox.sph/common.h"

namespace vox {
class Simulator_GUI_Base {
protected:
    SimulatorBase *m_simulatorBase;

public:
    Simulator_GUI_Base(SimulatorBase *simulatorBase) : m_simulatorBase(simulatorBase){};
    virtual ~Simulator_GUI_Base(){};

public:
    virtual void init(int argc, char **argv, const char *name) {}
    virtual void initSimulationParameterGUI() {}
    virtual void initParameterGUI() {}
    virtual void render() {}
    virtual void reset() {}
    virtual void update() {}
    virtual void cleanup() {}
    virtual void run() {}
    virtual void stop() {}
    virtual void addKeyFunc(char k, std::function<void()> const &func) {}

    SimulatorBase *getSimulatorBase() const { return m_simulatorBase; }
};
}  // namespace vox