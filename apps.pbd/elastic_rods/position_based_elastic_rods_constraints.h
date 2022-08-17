//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Eigen/Dense>

#include "apps.pbd/elastic_rods/position_based_elastic_rods_model.h"
#include "vox.pbd/constraints.h"

namespace vox {
class SimulationModel;

class GhostPointEdgeDistanceConstraint : public Constraint {
public:
    static int TYPE_ID;
    Real m_restLength{};

    GhostPointEdgeDistanceConstraint() : Constraint(3) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    bool initConstraint(PositionBasedElasticRodsModel &model,
                        unsigned int particle1,
                        unsigned int particle2,
                        unsigned int particle3);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class PerpendiculaBisectorConstraint : public Constraint {
public:
    static int TYPE_ID;

    PerpendiculaBisectorConstraint() : Constraint(3) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    bool initConstraint(SimulationModel &model, unsigned int particle1, unsigned int particle2, unsigned int particle3);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class DarbouxVectorConstraint : public Constraint {
public:
    static int TYPE_ID;
    Matrix3r m_dA;  // material frame A
    Matrix3r m_dB;  // material frame B

    DarbouxVectorConstraint() : Constraint(5) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(PositionBasedElasticRodsModel &model,
                                unsigned int particle1,
                                unsigned int particle2,
                                unsigned int particle3,
                                unsigned int particle4,
                                unsigned int particle5);

    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

}  // namespace vox