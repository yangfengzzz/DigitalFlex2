//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/binary_file_reader_writer.h"
#include "vox.base/common.h"
#include "vox.base/sph_kernels.h"
#include "vox.sph/rigid_body_object.h"

namespace vox {
class TimeStep;

/** \brief The boundary model stores the information required for boundary handling
 */
class BoundaryModel {
public:
    BoundaryModel();
    virtual ~BoundaryModel();

protected:
    RigidBodyObject *m_rigidBody{};
    std::vector<Vector3r> m_forcePerThread;
    std::vector<Vector3r> m_torquePerThread;

public:
    virtual void reset();

    virtual void performNeighborhoodSearchSort(){};

    virtual void saveState(BinaryFileWriter &binWriter){};
    virtual void loadState(BinaryFileReader &binReader){};

    RigidBodyObject *getRigidBodyObject() { return m_rigidBody; }

    FORCE_INLINE void addForce(const Vector3r &pos, const Vector3r &f) {
        if (m_rigidBody->isDynamic()) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            m_forcePerThread[tid] += f;
            m_torquePerThread[tid] += (pos - m_rigidBody->getPosition()).cross(f);
        }
    }

    FORCE_INLINE void getPointVelocity(const Vector3r &x, Vector3r &res) {
        if (m_rigidBody->isDynamic() || m_rigidBody->isAnimated())
            res = m_rigidBody->getAngularVelocity().cross(x - m_rigidBody->getPosition()) + m_rigidBody->getVelocity();
        else
            res.setZero();
    }

    void getForceAndTorque(Vector3r &force, Vector3r &torque);
    void clearForceAndTorque();
};
}  // namespace vox
