//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"

namespace vox {
/** \brief Base class for rigid body objects.
 */
class RigidBodyObject {
protected:
    bool m_isAnimated;

public:
    RigidBodyObject() { m_isAnimated = false; }
    virtual ~RigidBodyObject() = default;
    ;

    [[nodiscard]] virtual bool isDynamic() const = 0;
    [[nodiscard]] bool isAnimated() const { return m_isAnimated; }
    virtual void setIsAnimated(const bool b) { m_isAnimated = true; }

    [[nodiscard]] virtual Real const getMass() const = 0;
    [[nodiscard]] virtual Vector3r const &getPosition() const = 0;
    virtual void setPosition(const Vector3r &x) = 0;
    [[nodiscard]] virtual Vector3r getWorldSpacePosition() const = 0;
    [[nodiscard]] virtual Vector3r const &getVelocity() const = 0;
    virtual void setVelocity(const Vector3r &v) = 0;
    [[nodiscard]] virtual Quaternionr const &getRotation() const = 0;
    virtual void setRotation(const Quaternionr &q) = 0;
    [[nodiscard]] virtual Matrix3r getWorldSpaceRotation() const = 0;
    [[nodiscard]] virtual Vector3r const &getAngularVelocity() const = 0;
    virtual void setAngularVelocity(const Vector3r &v) = 0;
    virtual void addForce(const Vector3r &f) = 0;
    virtual void addTorque(const Vector3r &t) = 0;

    virtual void updateMeshTransformation() = 0;
    [[nodiscard]] virtual const std::vector<Vector3r> &getVertices() const = 0;
    [[nodiscard]] virtual const std::vector<Vector3r> &getVertexNormals() const = 0;
    [[nodiscard]] virtual const std::vector<unsigned int> &getFaces() const = 0;
};
}  // namespace vox