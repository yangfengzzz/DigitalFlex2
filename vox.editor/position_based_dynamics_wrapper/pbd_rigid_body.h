//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.pbd/rigid_body.h"
#include "vox.base/common.h"
#include "vox.sph/rigid_body_object.h"
#include "vox.sph/time_manager.h"

namespace vox {
class PBDRigidBody : public RigidBodyObject {
protected:
    RigidBody *m_rigidBody;

public:
    PBDRigidBody(RigidBody *rigidBody) : RigidBodyObject(), m_rigidBody(rigidBody) { m_isAnimated = false; }

    [[nodiscard]] bool isDynamic() const override { return m_rigidBody->getMass() != 0.0; }
    void setIsAnimated(const bool b) override {}

    [[nodiscard]] Real const getMass() const override { return m_rigidBody->getMass(); }
    [[nodiscard]] Vector3r const &getPosition() const override { return m_rigidBody->getPosition(); }
    [[nodiscard]] Vector3r const &getVelocity() const override { return m_rigidBody->getVelocity(); }
    void setVelocity(const Vector3r &v) override { m_rigidBody->setVelocity(v); }
    [[nodiscard]] Quaternionr const &getRotation() const override { return m_rigidBody->getRotation(); }
    [[nodiscard]] Vector3r const &getAngularVelocity() const override { return m_rigidBody->getAngularVelocity(); }
    void setAngularVelocity(const Vector3r &v) override { m_rigidBody->setAngularVelocity(v); }

    void setPosition(const Vector3r &x) override {
        m_rigidBody->setPosition(x);
        m_rigidBody->getGeometry().updateMeshTransformation(m_rigidBody->getPosition(),
                                                            m_rigidBody->getRotationMatrix());
    }

    void setRotation(const Quaternionr &q) override {
        m_rigidBody->setRotation(q);
        m_rigidBody->rotationUpdated();
        m_rigidBody->getGeometry().updateMeshTransformation(m_rigidBody->getPosition(),
                                                            m_rigidBody->getRotationMatrix());
    }

    void addForce(const Vector3r &f) override {
        const Real dt = TimeManager::getCurrent()->getTimeStepSize();
        m_rigidBody->getVelocity() += (1.0 / m_rigidBody->getMass()) * f * dt;
    }

    void addTorque(const Vector3r &t) override {
        const Real dt = TimeManager::getCurrent()->getTimeStepSize();
        m_rigidBody->getAngularVelocity() += m_rigidBody->getInertiaTensorInverseW() * t * dt;
    }

    // transformation local to:
    // p_world = R R_MAT^T (R_initial p_local + x_initial - x_MAT) + x
    [[nodiscard]] Vector3r getWorldSpacePosition() const override {
        const Matrix3r R2 = (m_rigidBody->getRotation() * m_rigidBody->getRotationMAT().inverse()).matrix();
        const Vector3r x = R2 * m_rigidBody->getPositionInitial_MAT() + m_rigidBody->getPosition();
        return x;
    }

    [[nodiscard]] Matrix3r getWorldSpaceRotation() const override {
        return (m_rigidBody->getRotation() * m_rigidBody->getRotationMAT().inverse() *
                m_rigidBody->getRotationInitial())
                .matrix();
    }

    void updateMeshTransformation() override {
        m_rigidBody->getGeometry().updateMeshTransformation(m_rigidBody->getPosition(),
                                                            m_rigidBody->getRotationMatrix());
    }

    [[nodiscard]] const std::vector<Vector3r> &getVertices() const override {
        return m_rigidBody->getGeometry().getVertexData().getVertices();
    };
    [[nodiscard]] const std::vector<Vector3r> &getVertexNormals() const override {
        return m_rigidBody->getGeometry().getMesh().getVertexNormals();
    };
    [[nodiscard]] const std::vector<unsigned int> &getFaces() const override {
        return m_rigidBody->getGeometry().getMesh().getFaces();
    };
};
}  // namespace vox