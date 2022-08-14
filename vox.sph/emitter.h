//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"

namespace vox {
class Emitter {
public:
    Emitter(FluidModel *model,
            unsigned int width,
            unsigned int height,
            Vector3r pos,
            Matrix3r rotation,
            Real velocity,
            unsigned int type = 0);
    virtual ~Emitter();

protected:
    FluidModel *m_model;
    unsigned int m_width;
    unsigned int m_height;
    Vector3r m_x;
    Matrix3r m_rotation;
    Real m_velocity;
    unsigned int m_type;
    Real m_nextEmitTime;
    Real m_emitStartTime;
    Real m_emitEndTime;
    unsigned int m_emitCounter;
    unsigned int m_objectId{};

    static FORCE_INLINE bool inBox(const Vector3r &x,
                                   const Vector3r &xBox,
                                   const Matrix3r &rotBox,
                                   const Vector3r &scaleBox) {
        const Vector3r xlocal = rotBox.transpose() * (x - xBox);
        // for a box shape, m_scale stores the half-size of the box
        // inside box if closer than half-size on all axes
        return (xlocal.array().abs() < scaleBox.array()).all();
    }

    static FORCE_INLINE bool inCylinder(
            const Vector3r &x, const Vector3r &xCyl, const Matrix3r &rotCyl, const Real h, const Real r2) {
        const Vector3r xlocal = rotCyl.transpose() * (x - xCyl);
        // inside cylinder if distance to x-axis is less than r
        // and projection on x-axis is between 0 and h
        const Real proj = xlocal.x();
        const Real d2 = Vector2r(xlocal.y(), xlocal.z()).squaredNorm();
        const Real hHalf = static_cast<Real>(0.5) * h;
        return (proj > -hHalf) && (proj < hHalf) && (d2 < r2);
    }

public:
    void emitParticles(std::vector<unsigned int> &reusedParticles,
                       unsigned int &indexReuse,
                       unsigned int &numEmittedParticles);
    void emitParticlesCircle(std::vector<unsigned int> &reusedParticles,
                             unsigned int &indexReuse,
                             unsigned int &numEmittedParticles);
    [[nodiscard]] Real getNextEmitTime() const { return m_nextEmitTime; }
    void setNextEmitTime(Real val) { m_nextEmitTime = val; }
    void setEmitStartTime(Real val) {
        m_emitStartTime = val;
        setNextEmitTime(val);
    }
    void setEmitEndTime(Real val) { m_emitEndTime = val; }
    static Vector3r getSize(Real width, Real height, int type);

    void step(std::vector<unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles);
    virtual void reset();

    void saveState(BinaryFileWriter &binWriter) const;
    void loadState(BinaryFileReader &binReader);

    [[nodiscard]] const Vector3r &getPosition() const { return m_x; }
    void setPosition(const Vector3r &x) { m_x = x; }
    [[nodiscard]] const Matrix3r &getRotation() const { return m_rotation; }
    void setRotation(const Matrix3r &r) { m_rotation = r; }
    [[nodiscard]] Real getVelocity() const { return m_velocity; }
    void setVelocity(const Real v) { m_velocity = v; }
    [[nodiscard]] unsigned int getObjectId() const { return m_objectId; }
    void setObjectId(const unsigned int v) { m_objectId = v; }
};
}  // namespace vox