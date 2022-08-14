//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/drag/drag_force_macklin2014.h"

using namespace vox;

DragForce_Macklin2014::DragForce_Macklin2014(FluidModel *model) : DragBase(model) {}

DragForce_Macklin2014::~DragForce_Macklin2014() = default;

void DragForce_Macklin2014::step() {
    const Real density0 = m_model->getDensity0();

    const unsigned int numParticles = m_model->numActiveParticles();

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            Vector3r &ai = m_model->getAcceleration(i);
            const Vector3r &vi = m_model->getVelocity(i);
            ai -= m_dragCoefficient * static_cast<Real>(1.0) / m_model->getMass(i) * vi *
                  (1.0 - m_model->getDensity(i) / density0);
        }
    }
}

void DragForce_Macklin2014::reset() {}
