//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/rigid_body_object.h"
#include "vox.sph/time_manager.h"
#include "vox.sph/triangle_mesh.h"

namespace vox {
/** \brief This class stores the information of a static rigid body which
 * is not part of a rigid body simulation.
 */
class StaticRigidBody : public RigidBodyObject {
protected:
    Vector3r m_x0;
    Vector3r m_x;
    Quaternionr m_q;
    Quaternionr m_q0;
    Vector3r m_velocity;
    Vector3r m_angularVelocity;
    TriangleMesh m_geometry;

public:
    StaticRigidBody() {
        m_isAnimated = false;
        m_velocity = Vector3r::Zero();
        m_angularVelocity = Vector3r::Zero();
    }

    [[nodiscard]] bool isDynamic() const override { return false; }

    [[nodiscard]] Real const getMass() const override { return 0.0; }
    [[nodiscard]] Vector3r const& getPosition() const override { return m_x; }
    void setPosition(const Vector3r& x) override { m_x = x; }
    [[nodiscard]] Vector3r const& getPosition0() const { return m_x0; }
    void setPosition0(const Vector3r& x) { m_x0 = x; }
    [[nodiscard]] Vector3r getWorldSpacePosition() const override { return m_x; }
    [[nodiscard]] Vector3r const& getVelocity() const override { return m_velocity; }
    void setVelocity(const Vector3r& v) override {
        if (m_isAnimated) m_velocity = v;
    }
    [[nodiscard]] Quaternionr const& getRotation() const override { return m_q; }
    void setRotation(const Quaternionr& q) override { m_q = q; }
    [[nodiscard]] Quaternionr const& getRotation0() const { return m_q0; }
    void setRotation0(const Quaternionr& q) { m_q0 = q; }
    [[nodiscard]] Matrix3r getWorldSpaceRotation() const override { return m_q.toRotationMatrix(); }
    [[nodiscard]] Vector3r const& getAngularVelocity() const override { return m_angularVelocity; }
    void setAngularVelocity(const Vector3r& v) override {
        if (m_isAnimated) m_angularVelocity = v;
    }
    void addForce(const Vector3r& f) override {}
    void addTorque(const Vector3r& t) override {}
    void animate() {
        const Real dt = TimeManager::getCurrent()->getTimeStepSize();
        m_x += m_velocity * dt;
        Quaternionr angVelQ(0.0, m_angularVelocity[0], m_angularVelocity[1], m_angularVelocity[2]);
        m_q.coeffs() += dt * 0.5 * (angVelQ * m_q).coeffs();
        m_q.normalize();
        updateMeshTransformation();
    }

    [[nodiscard]] const std::vector<Vector3r>& getVertices() const override { return m_geometry.getVertices(); }
    [[nodiscard]] const std::vector<Vector3r>& getVertexNormals() const override {
        return m_geometry.getVertexNormals();
    }
    [[nodiscard]] const std::vector<unsigned int>& getFaces() const override { return m_geometry.getFaces(); }

    void setWorldSpacePosition(const Vector3r& x) { m_x = x; }
    void setWorldSpaceRotation(const Matrix3r& r) { m_q = Quaternionr(r); }
    TriangleMesh& getGeometry() { return m_geometry; }

    void updateMeshTransformation() override {
        m_geometry.updateMeshTransformation(m_x, m_q.toRotationMatrix());
        m_geometry.updateNormals();
        m_geometry.updateVertexNormals();
    }

    void reset() {
        m_x = m_x0;
        m_q = m_q0;
        updateMeshTransformation();
    }
};
}  // namespace vox