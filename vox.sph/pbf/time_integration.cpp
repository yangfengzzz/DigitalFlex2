//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/pbf/time_integration.h"

using namespace vox;

// ----------------------------------------------------------------------------------------------
void TimeIntegration::semiImplicitEuler(
        const Real h, const Real mass, Vector3r &position, Vector3r &velocity, const Vector3r &acceleration) {
    if (mass != 0.0) {
        velocity += acceleration * h;
        position += velocity * h;
    }
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::velocityUpdateFirstOrder(
        const Real h, const Real mass, const Vector3r &position, const Vector3r &oldPosition, Vector3r &velocity) {
    if (mass != 0.0) velocity = (1.0 / h) * (position - oldPosition);
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::velocityUpdateSecondOrder(const Real h,
                                                const Real mass,
                                                const Vector3r &position,
                                                const Vector3r &oldPosition,
                                                const Vector3r &positionOfLastStep,
                                                Vector3r &velocity) {
    if (mass != 0.0) velocity = (1.0 / h) * (1.5 * position - 2.0 * oldPosition + 0.5 * positionOfLastStep);
}
