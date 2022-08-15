//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/boundary_model_bender2019.h"

#include <iostream>

#include "vox.base/logging.h"
#include "vox.base/compact_search/compact_search.h"
#include "vox.sph/simulation.h"
#include "vox.sph/time_step.h"

using namespace vox;

BoundaryModel_Bender2019::BoundaryModel_Bender2019() : m_boundaryVolume(), m_boundaryXj() {
    m_map = nullptr;
    m_maxDist = 0.0;
    m_maxVel = 0.0;
}

BoundaryModel_Bender2019::~BoundaryModel_Bender2019() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nModels = sim->numberOfFluidModels();
    for (unsigned int i = 0; i < nModels; i++) {
        m_boundaryVolume[i].clear();
        m_boundaryXj[i].clear();
    }
    m_boundaryVolume.clear();
    m_boundaryXj.clear();

    delete m_map;
}

void BoundaryModel_Bender2019::initModel(RigidBodyObject *rbo) {
    Simulation *sim = Simulation::getCurrent();

    const unsigned int nModels = sim->numberOfFluidModels();
    m_boundaryVolume.resize(nModels);
    m_boundaryXj.resize(nModels);
    for (unsigned int i = 0; i < nModels; i++) {
        FluidModel *fm = sim->getFluidModel(i);
        m_boundaryVolume[i].resize(fm->numParticles(), 0.0);
        m_boundaryXj[i].resize(fm->numParticles(), Vector3r::Zero());
    }

    if (rbo->isDynamic()) {
#ifdef _OPENMP
        const int maxThreads = omp_get_max_threads();
#else
        const int maxThreads = 1;
#endif
        m_forcePerThread.resize(maxThreads, Vector3r::Zero());
        m_torquePerThread.resize(maxThreads, Vector3r::Zero());
    }

    m_rigidBody = rbo;
}

void BoundaryModel_Bender2019::reset() {
    BoundaryModel::reset();

    m_maxDist = 0.0;
    m_maxVel = 0.0;
}
