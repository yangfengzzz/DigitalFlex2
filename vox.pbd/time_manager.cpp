//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.pbd/time_manager.h"

using namespace vox;

TimeManager *TimeManager::current = nullptr;

TimeManager::TimeManager() {
    time = 0;
    h = static_cast<Real>(0.005);
}

TimeManager::~TimeManager() { current = 0; }

TimeManager *TimeManager::getCurrent() {
    if (current == 0) {
        current = new TimeManager();
    }
    return current;
}

void TimeManager::setCurrent(TimeManager *tm) { current = tm; }

bool TimeManager::hasCurrent() { return (current != 0); }

Real TimeManager::getTime() { return time; }

void TimeManager::setTime(Real t) { time = t; }

Real TimeManager::getTimeStepSize() { return h; }

void TimeManager::setTimeStepSize(Real tss) { h = tss; }
