//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"
#include "vox.sph/gui//simulator_gui_base.h"

namespace vox {
class Simulator_GUI_imgui : public Simulator_GUI_Base {
public:
    Simulator_GUI_imgui(SimulatorBase *base);
    virtual ~Simulator_GUI_imgui();

protected:
    unsigned int m_currentFluidModel;
    Vector3r m_oldMousePos;
    std::vector<std::string> m_colorFieldNames;
    std::vector<std::vector<unsigned int>> m_selectedParticles;

    std::vector<std::vector<unsigned int>> &getSelectedParticles() { return m_selectedParticles; }
    void initImgui();
    void initImguiParameters();
    void renderBoundary();
    void particleInfo();

    static void selection(const Vector2i &start, const Vector2i &end, void *clientData);
    static void mouseMove(int x, int y, void *clientData);

    void switchPause();
    static void switchDrawMode();

    void destroy();

public:
    virtual void init(int argc, char **argv, const char *name);
    virtual void initParameterGUI();
    virtual void initSimulationParameterGUI();
    virtual void render();
    virtual void reset();
    virtual void update();
    virtual void cleanup();
    virtual void run();
    virtual void stop();
    virtual void addKeyFunc(char k, std::function<void()> const &func);

    void createSimulationParameterGUI();
};
}  // namespace vox