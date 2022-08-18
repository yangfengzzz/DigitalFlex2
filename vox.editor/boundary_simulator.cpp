//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.editor/boundary_simulator.h"

#include "vox.base/time_manager.h"
#include "vox.sph/boundary_model.h"
#include "vox.sph/simulation.h"

using namespace vox;

void BoundarySimulator::updateBoundaryForces() {
    Real h = TimeManager::getCurrent()->getTimeStepSize();
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nObjects = sim->numberOfBoundaryModels();
    for (unsigned int i = 0; i < nObjects; i++) {
        BoundaryModel *bm = sim->getBoundaryModel(i);
        RigidBodyObject *rbo = bm->getRigidBodyObject();
        if (rbo->isDynamic()) {
            Vector3r force, torque;
            bm->getForceAndTorque(force, torque);
            rbo->addForce(force);
            rbo->addTorque(torque);
            bm->clearForceAndTorque();
        }
    }
}
