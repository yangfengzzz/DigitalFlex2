//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/eigen.h>
#include <pybind11/stl_bind.h>

#include "py.pbd/common.h"

namespace py = pybind11;

#define CONSTRAINT(name, base)                     \
    py::class_<vox::name, vox::base>(m_sub, #name) \
            .def(py::init<>())                     \
            .def_readwrite_static("TYPE_ID", &vox::name::TYPE_ID)

#define CONSTRAINT_JOINTINFO(name, base) CONSTRAINT(name, base).def_readwrite("jointInfo", &vox::name::m_jointInfo)

void ConstraintsModule(py::module m_sub) {
    py::class_<vox::Constraint>(m_sub, "Constraint")
            .def_readwrite("bodies", &vox::Constraint::m_bodies)
            .def("getTypeId", &vox::Constraint::getTypeId)
            .def("initConstraintBeforeProjection", &vox::Constraint::initConstraintBeforeProjection)
            .def("updateConstraint", &vox::Constraint::updateConstraint)
            .def("solvePositionConstraint", &vox::Constraint::solvePositionConstraint)
            .def("solveVelocityConstraint", &vox::Constraint::solveVelocityConstraint);

    CONSTRAINT_JOINTINFO(BallJoint, Constraint);
    CONSTRAINT_JOINTINFO(BallOnLineJoint, Constraint);
    CONSTRAINT_JOINTINFO(HingeJoint, Constraint);
    CONSTRAINT_JOINTINFO(UniversalJoint, Constraint);
    CONSTRAINT_JOINTINFO(SliderJoint, Constraint);

    py::class_<vox::MotorJoint, vox::Constraint>(m_sub, "MotorJoint")
            .def("getTarget", &vox::MotorJoint::getTarget)
            .def("setTarget", &vox::MotorJoint::setTarget)
            .def("getTargetSequence", &vox::MotorJoint::getTargetSequence)
            .def("setTargetSequence", &vox::MotorJoint::setTargetSequence)
            .def("getRepeatSequence", &vox::MotorJoint::getRepeatSequence)
            .def("setRepeatSequence", &vox::MotorJoint::setRepeatSequence);

    CONSTRAINT_JOINTINFO(TargetPositionMotorSliderJoint, MotorJoint);
    CONSTRAINT_JOINTINFO(TargetVelocityMotorSliderJoint, MotorJoint);
    CONSTRAINT_JOINTINFO(TargetAngleMotorHingeJoint, MotorJoint);
    CONSTRAINT_JOINTINFO(TargetVelocityMotorHingeJoint, MotorJoint);
    CONSTRAINT_JOINTINFO(DamperJoint, Constraint)
            .def_readwrite("stiffness", &vox::DamperJoint::m_stiffness)
            .def_readwrite("lambda", &vox::DamperJoint::m_lambda);
    CONSTRAINT_JOINTINFO(RigidBodyParticleBallJoint, Constraint);
    CONSTRAINT_JOINTINFO(RigidBodySpring, Constraint)
            .def_readwrite("restLength", &vox::RigidBodySpring::m_restLength)
            .def_readwrite("stiffness", &vox::RigidBodySpring::m_stiffness)
            .def_readwrite("lambda", &vox::RigidBodySpring::m_lambda);
    CONSTRAINT_JOINTINFO(DistanceJoint, Constraint).def_readwrite("restLength", &vox::DistanceJoint::m_restLength);

    CONSTRAINT(DistanceConstraint, Constraint)
            .def_readwrite("stiffness", &vox::DistanceConstraint::m_stiffness)
            .def_readwrite("restLength", &vox::DistanceConstraint::m_restLength);
    CONSTRAINT(DistanceConstraint_XPBD, Constraint)
            .def_readwrite("stiffness", &vox::DistanceConstraint_XPBD::m_stiffness)
            .def_readwrite("restLength", &vox::DistanceConstraint_XPBD::m_restLength)
            .def_readwrite("lambda", &vox::DistanceConstraint_XPBD::m_lambda);
    CONSTRAINT(DihedralConstraint, Constraint)
            .def_readwrite("stiffness", &vox::DihedralConstraint::m_stiffness)
            .def_readwrite("restAngle", &vox::DihedralConstraint::m_restAngle);
    CONSTRAINT(IsometricBendingConstraint, Constraint)
            .def_readwrite("stiffness", &vox::IsometricBendingConstraint::m_stiffness)
            .def_readwrite("Q", &vox::IsometricBendingConstraint::m_Q);
    CONSTRAINT(IsometricBendingConstraint_XPBD, Constraint)
            .def_readwrite("stiffness", &vox::IsometricBendingConstraint_XPBD::m_stiffness)
            .def_readwrite("Q", &vox::IsometricBendingConstraint_XPBD::m_Q)
            .def_readwrite("lambda", &vox::IsometricBendingConstraint_XPBD::m_lambda);
    CONSTRAINT(FEMTriangleConstraint, Constraint)
            .def_readwrite("xxStiffness", &vox::FEMTriangleConstraint::m_xxStiffness)
            .def_readwrite("yyStiffness", &vox::FEMTriangleConstraint::m_yyStiffness)
            .def_readwrite("xyStiffness", &vox::FEMTriangleConstraint::m_xyStiffness)
            .def_readwrite("xyPoissonRatio", &vox::FEMTriangleConstraint::m_xyPoissonRatio)
            .def_readwrite("yxPoissonRatio", &vox::FEMTriangleConstraint::m_yxPoissonRatio)
            .def_readwrite("area", &vox::FEMTriangleConstraint::m_area)
            .def_readwrite("invRestMat", &vox::FEMTriangleConstraint::m_invRestMat);
    CONSTRAINT(StrainTriangleConstraint, Constraint)
            .def_readwrite("xxStiffness", &vox::StrainTriangleConstraint::m_xxStiffness)
            .def_readwrite("yyStiffness", &vox::StrainTriangleConstraint::m_yyStiffness)
            .def_readwrite("xyStiffness", &vox::StrainTriangleConstraint::m_xyStiffness)
            .def_readwrite("normalizeStretch", &vox::StrainTriangleConstraint::m_normalizeStretch)
            .def_readwrite("normalizeShear", &vox::StrainTriangleConstraint::m_normalizeShear)
            .def_readwrite("invRestMat", &vox::StrainTriangleConstraint::m_invRestMat);
    CONSTRAINT(VolumeConstraint, Constraint)
            .def_readwrite("stiffness", &vox::VolumeConstraint::m_stiffness)
            .def_readwrite("restVolume", &vox::VolumeConstraint::m_restVolume);
    CONSTRAINT(VolumeConstraint_XPBD, Constraint)
            .def_readwrite("stiffness", &vox::VolumeConstraint_XPBD::m_stiffness)
            .def_readwrite("restVolume", &vox::VolumeConstraint_XPBD::m_restVolume)
            .def_readwrite("lambda", &vox::VolumeConstraint_XPBD::m_lambda);
    CONSTRAINT(FEMTetConstraint, Constraint)
            .def_readwrite("stiffness", &vox::FEMTetConstraint::m_stiffness)
            .def_readwrite("poissonRatio", &vox::FEMTetConstraint::m_poissonRatio)
            .def_readwrite("volume", &vox::FEMTetConstraint::m_volume)
            .def_readwrite("invRestMat", &vox::FEMTetConstraint::m_invRestMat);
    CONSTRAINT(StrainTetConstraint, Constraint)
            .def_readwrite("stretchStiffness", &vox::StrainTetConstraint::m_stretchStiffness)
            .def_readwrite("shearStiffness", &vox::StrainTetConstraint::m_shearStiffness)
            .def_readwrite("normalizeStretch", &vox::StrainTetConstraint::m_normalizeStretch)
            .def_readwrite("normalizeShear", &vox::StrainTetConstraint::m_normalizeShear)
            .def_readwrite("invRestMat", &vox::StrainTetConstraint::m_invRestMat);
    py::class_<vox::ShapeMatchingConstraint, vox::Constraint>(m_sub, "ShapeMatchingConstraint")
            .def(py::init<const unsigned int>())
            .def_readwrite("stiffness", &vox::ShapeMatchingConstraint::m_stiffness)
            .def_readwrite_static("TYPE_ID", &vox::ShapeMatchingConstraint::TYPE_ID)
            .def_readwrite("restCm", &vox::ShapeMatchingConstraint::m_restCm);
    CONSTRAINT(StretchShearConstraint, Constraint)
            .def_readwrite("stretchingStiffness", &vox::StretchShearConstraint::m_stretchingStiffness)
            .def_readwrite("shearingStiffness1", &vox::StretchShearConstraint::m_shearingStiffness1)
            .def_readwrite("shearingStiffness2", &vox::StretchShearConstraint::m_shearingStiffness2)
            .def_readwrite("restLength", &vox::StretchShearConstraint::m_restLength);
    CONSTRAINT(BendTwistConstraint, Constraint)
            .def_readwrite("twistingStiffness", &vox::BendTwistConstraint::m_twistingStiffness)
            .def_readwrite("bendingStiffness1", &vox::BendTwistConstraint::m_bendingStiffness1)
            .def_readwrite("bendingStiffness2", &vox::BendTwistConstraint::m_bendingStiffness2)
            .def_readwrite("restDarbouxVector", &vox::BendTwistConstraint::m_restDarbouxVector);
    CONSTRAINT(StretchBendingTwistingConstraint, Constraint)
            .def_readwrite("averageRadius", &vox::StretchBendingTwistingConstraint::m_averageRadius)
            .def_readwrite("averageSegmentLength", &vox::StretchBendingTwistingConstraint::m_averageSegmentLength)
            .def_readwrite("restDarbouxVector", &vox::StretchBendingTwistingConstraint::m_restDarbouxVector)
            .def_readwrite("stiffnessCoefficientK", &vox::StretchBendingTwistingConstraint::m_stiffnessCoefficientK)
            .def_readwrite("stretchCompliance", &vox::StretchBendingTwistingConstraint::m_stretchCompliance)
            .def_readwrite("bendingAndTorsionCompliance",
                           &vox::StretchBendingTwistingConstraint::m_bendingAndTorsionCompliance)
            .def_readwrite("lambdaSum", &vox::StretchBendingTwistingConstraint::m_lambdaSum);
    CONSTRAINT(DirectPositionBasedSolverForStiffRodsConstraint, Constraint)
            .def("initConstraint", &vox::DirectPositionBasedSolverForStiffRodsConstraint::initConstraint);

    py::class_<vox::RigidBodyContactConstraint>(m_sub, "RigidBodyContactConstraint")
            .def_readwrite("bodies", &vox::RigidBodyContactConstraint::m_bodies)
            .def_readwrite("stiffness", &vox::RigidBodyContactConstraint::m_stiffness)
            .def_readwrite("frictionCoeff", &vox::RigidBodyContactConstraint::m_frictionCoeff)
            .def_readwrite("sum_impulses", &vox::RigidBodyContactConstraint::m_sum_impulses)
            .def_readwrite("constraintInfo", &vox::RigidBodyContactConstraint::m_constraintInfo)
            .def("getTypeId", &vox::RigidBodyContactConstraint::getTypeId)
            .def("initConstraint", &vox::RigidBodyContactConstraint::initConstraint)
            .def("solveVelocityConstraint", &vox::RigidBodyContactConstraint::solveVelocityConstraint);

    py::class_<vox::ParticleRigidBodyContactConstraint>(m_sub, "ParticleRigidBodyContactConstraint")
            .def_readwrite("bodies", &vox::ParticleRigidBodyContactConstraint::m_bodies)
            .def_readwrite("stiffness", &vox::ParticleRigidBodyContactConstraint::m_stiffness)
            .def_readwrite("frictionCoeff", &vox::ParticleRigidBodyContactConstraint::m_frictionCoeff)
            .def_readwrite("sum_impulses", &vox::ParticleRigidBodyContactConstraint::m_sum_impulses)
            .def_readwrite("constraintInfo", &vox::ParticleRigidBodyContactConstraint::m_constraintInfo)
            .def("initConstraint", &vox::ParticleRigidBodyContactConstraint::initConstraint)
            .def("solveVelocityConstraint", &vox::ParticleRigidBodyContactConstraint::solveVelocityConstraint);

    py::class_<vox::ParticleTetContactConstraint>(m_sub, "ParticleTetContactConstraint")
            .def_readwrite("bodies", &vox::ParticleTetContactConstraint::m_bodies)
            .def_readwrite("solidIndex", &vox::ParticleTetContactConstraint::m_solidIndex)
            .def_readwrite("tetIndex", &vox::ParticleTetContactConstraint::m_tetIndex)
            .def_readwrite("frictionCoeff", &vox::ParticleTetContactConstraint::m_frictionCoeff)
            .def_readwrite("bary", &vox::ParticleTetContactConstraint::m_bary)
            .def_readwrite("lambda", &vox::ParticleTetContactConstraint::m_lambda)
            .def_readwrite("constraintInfo", &vox::ParticleTetContactConstraint::m_constraintInfo)
            .def_readwrite("x", &vox::ParticleTetContactConstraint::m_x)
            .def_readwrite("v", &vox::ParticleTetContactConstraint::m_v)
            .def("initConstraint", &vox::ParticleTetContactConstraint::initConstraint)
            .def("solvePositionConstraint", &vox::ParticleTetContactConstraint::solvePositionConstraint)
            .def("solveVelocityConstraint", &vox::ParticleTetContactConstraint::solveVelocityConstraint);
}