//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.pbd/simulation.h"

#include "vox.pbd/time_manager.h"
#include "vox.pbd/time_step_controller.h"
#include "vox.pbd/timing.h"

using namespace vox;
using namespace std;

simulation *simulation::current = nullptr;
int simulation::GRAVITATION = -1;

simulation::simulation() {
    m_gravitation = Vector3r(0.0, -9.81, 0.0);

    m_timeStep = nullptr;
}

simulation::~simulation() {
    delete m_timeStep;
    delete TimeManager::getCurrent();

    current = nullptr;
}

simulation *simulation::getCurrent() {
    if (current == nullptr) {
        current = new simulation();
        current->init();
    }
    return current;
}

void simulation::setCurrent(simulation *tm) { current = tm; }

bool simulation::hasCurrent() { return (current != nullptr); }

void simulation::init() {
    initParameters();

    m_timeStep = new TimeStepController();
    m_timeStep->init();
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));
}

void simulation::initParameters() {
    ParameterObject::initParameters();

    GRAVITATION = createVectorParameter("gravitation", "Gravitation", 3u, m_gravitation.data());
    setGroup(GRAVITATION, "simulation");
    setDescription(GRAVITATION, "Vector to define the gravitational acceleration.");
}

void simulation::reset() {
    m_model->reset();
    if (m_timeStep) m_timeStep->reset();

    TimeManager::getCurrent()->setTime(static_cast<Real>(0.0));
}
