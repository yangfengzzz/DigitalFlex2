//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Eigen/Dense>

#include "vox.pbd/constraints.h"

namespace vox {
class SimulationModel;

class GenericDistanceConstraint : public Constraint {
public:
    static int TYPE_ID;
    Real m_restLength{};
    Real m_stiffness{};

    static void constraintFct(unsigned int numberOfParticles,
                              const Real invMass[],
                              const Vector3r x[],
                              void *userData,
                              Eigen::Matrix<Real, 1, 1> &constraintValue);

    static void gradientFct(unsigned int i,
                            unsigned int numberOfParticles,
                            const Real invMass[],
                            const Vector3r x[],
                            void *userData,
                            Eigen::Matrix<Real, 1, 3> &jacobian);

    GenericDistanceConstraint() : Constraint(2) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(SimulationModel &model,
                                unsigned int particle1,
                                unsigned int particle2,
                                Real stiffness);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class GenericIsometricBendingConstraint : public Constraint {
public:
    static int TYPE_ID;
    Matrix4r m_Q;
    Real m_stiffness{};

    static void constraintFct(unsigned int numberOfParticles,
                              const Real invMass[],
                              const Vector3r x[],
                              void *userData,
                              Eigen::Matrix<Real, 1, 1> &constraintValue);

    GenericIsometricBendingConstraint() : Constraint(4) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(SimulationModel &model,
                                unsigned int particle1,
                                unsigned int particle2,
                                unsigned int particle3,
                                unsigned int particle4,
                                Real stiffness);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class GenericHingeJoint : public Constraint {
public:
    static int TYPE_ID;
    Eigen::Matrix<Real, 3, 12> m_jointInfo;

    static void constraintFct(unsigned int numberOfRigidBodies,
                              const Real invMass[],              // inverse mass is zero if body is static
                              const Vector3r x[],                // positions of bodies
                              const Matrix3r inertiaInverseW[],  // inverse inertia tensor (world space) of bodies
                              const Quaternionr q[],
                              void *userData,
                              Eigen::Matrix<Real, 5, 1> &constraintValue);

    static void gradientFct(unsigned int i,
                            unsigned int numberOfRigidBodies,
                            const Real invMass[],
                            const Vector3r x[],
                            const Matrix3r inertiaInverseW[],
                            const Quaternionr q[],
                            void *userData,
                            Eigen::Matrix<Real, 5, 6> &jacobian);

    GenericHingeJoint() : Constraint(2) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(SimulationModel &model,
                                unsigned int rbIndex1,
                                unsigned int rbIndex2,
                                const Vector3r &pos,
                                const Vector3r &axis);
    bool updateConstraint(SimulationModel &model) override;
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class GenericBallJoint : public Constraint {
public:
    static int TYPE_ID;
    Eigen::Matrix<Real, 3, 2> m_jointInfo;

    static void constraintFct(unsigned int numberOfRigidBodies,
                              const Real invMass[],              // inverse mass is zero if body is static
                              const Vector3r x[],                // positions of bodies
                              const Matrix3r inertiaInverseW[],  // inverse inertia tensor (world space) of bodies
                              const Quaternionr q[],
                              void *userData,
                              Eigen::Matrix<Real, 3, 1> &constraintValue);

    GenericBallJoint() : Constraint(2) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(SimulationModel &model,
                                unsigned int rbIndex1,
                                unsigned int rbIndex2,
                                const Vector3r &pos);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};

class GenericSliderJoint : public Constraint {
public:
    static int TYPE_ID;
    Eigen::Matrix<Real, 3, 7> m_jointInfo;

    static void constraintFct(unsigned int numberOfRigidBodies,
                              const Real invMass[],              // inverse mass is zero if body is static
                              const Vector3r x[],                // positions of bodies
                              const Matrix3r inertiaInverseW[],  // inverse inertia tensor (world space) of bodies
                              const Quaternionr q[],
                              void *userData,
                              Eigen::Matrix<Real, 5, 1> &constraintValue);

    GenericSliderJoint() : Constraint(2) {}
    [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }

    virtual bool initConstraint(SimulationModel &model,
                                unsigned int rbIndex1,
                                unsigned int rbIndex2,
                                const Vector3r &pos,
                                const Vector3r &axis);
    bool solvePositionConstraint(SimulationModel &model, unsigned int iter) override;
};
}  // namespace vox