//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.pbd/particle_data.h"
#include "vox.pbd/simulation_model.h"

namespace vox {
class Constraint;

class PositionBasedElasticRodsModel : public SimulationModel {
public:
    PositionBasedElasticRodsModel();
    ~PositionBasedElasticRodsModel() override;

protected:
    ParticleData m_ghostParticles;
    Vector3r m_restDarbouxVector;
    Vector3r m_stiffness;

public:
    virtual void reset();
    virtual void cleanup();

    ParticleData &getGhostParticles();
    void addElasticRodModel(unsigned int nPoints, Vector3r *points);

    bool addPerpendiculaBisectorConstraint(unsigned int p0, unsigned int p1, unsigned int p2);
    bool addGhostPointEdgeDistanceConstraint(unsigned int pA, unsigned int pB, unsigned int pG);
    bool addDarbouxVectorConstraint(
            unsigned int pA, unsigned int pB, unsigned int pC, unsigned int pD, unsigned int pE);

    void setRestDarbouxVector(const Vector3r &val) { m_restDarbouxVector = val; }
    Vector3r &getRestDarbouxVector() { return m_restDarbouxVector; }
    void setBendingAndTwistingStiffness(const Vector3r &val) { m_stiffness = val; }
    Vector3r &getBendingAndTwistingStiffness() { return m_stiffness; }
};
}  // namespace vox