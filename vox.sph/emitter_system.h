//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"
#include "vox.sph/emitter.h"

namespace vox {
class TimeStep;
class FluidModel;

class EmitterSystem {
public:
    explicit EmitterSystem(FluidModel *model);
    virtual ~EmitterSystem();

protected:
    FluidModel *m_model;
    static const unsigned int m_maxParticlesToReusePerStep = 50000;
    bool m_reuseParticles;
    Vector3r m_boxMin;
    Vector3r m_boxMax;
    unsigned int m_numberOfEmittedParticles;
    unsigned int m_numReusedParticles;
    std::vector<unsigned int> m_reusedParticles;
    std::vector<Emitter *> m_emitters;

    void reuseParticles();

public:
    void enableReuseParticles(const Vector3r &boxMin = Vector3r(-1, -1, -1),
                              const Vector3r &boxMax = Vector3r(1, 1, 1));
    void disableReuseParticles();
    void addEmitter(unsigned int width,
                    unsigned int height,
                    const Vector3r &pos,
                    const Matrix3r &rotation,
                    Real velocity,
                    unsigned int type);
    [[nodiscard]] unsigned int numEmitters() const { return static_cast<unsigned int>(m_emitters.size()); }
    std::vector<Emitter *> &getEmitters() { return m_emitters; }

    [[nodiscard]] unsigned int numReusedParticles() const { return m_numReusedParticles; }
    [[nodiscard]] unsigned int numEmittedParticles() const { return m_numberOfEmittedParticles; }

    void step();
    void reset();

    void saveState(BinaryFileWriter &binWriter);
    void loadState(BinaryFileReader &binReader);
};
}  // namespace vox