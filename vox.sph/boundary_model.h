//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/binary_file_reader_writer.h"
#include "vox.sph/common.h"
#include "vox.sph/rigid_body_object.h"
#include "vox.sph/sph_kernels.h"

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

#ifdef USE_AVX
    FORCE_INLINE void addForce(const Vector3f8 &pos, const Vector3f8 &f, const unsigned int count) {
        if (m_rigidBody->isDynamic()) {
            float fx[8];
            float fy[8];
            float fz[8];
            f.x().store(fx);
            f.y().store(fy);
            f.z().store(fz);
            float px[8];
            float py[8];
            float pz[8];
            pos.x().store(px);
            pos.y().store(py);
            pos.z().store(pz);
            for (unsigned int l = 0; l < count; l++) {
                addForce(Vector3r(px[l], py[l], pz[l]), Vector3r(fx[l], fy[l], fz[l]));
            }
        }
    }
#endif

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