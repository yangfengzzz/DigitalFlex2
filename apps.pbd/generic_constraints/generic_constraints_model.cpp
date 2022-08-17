//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "apps.pbd/generic_constraints/generic_constraints_model.h"

#include "apps.pbd/generic_constraints/generic_constraints.h"

using namespace vox;

GenericConstraintsModel::GenericConstraintsModel() : SimulationModel() {}

GenericConstraintsModel::~GenericConstraintsModel() = default;

bool GenericConstraintsModel::addGenericDistanceConstraint(const unsigned int particle1,
                                                           const unsigned int particle2,
                                                           const Real stiffness) {
    auto *c = new GenericDistanceConstraint();
    const bool res = c->initConstraint(*this, particle1, particle2, stiffness);
    if (res) m_constraints.push_back(c);
    return res;
}

bool GenericConstraintsModel::addGenericIsometricBendingConstraint(const unsigned int particle1,
                                                                   const unsigned int particle2,
                                                                   const unsigned int particle3,
                                                                   const unsigned int particle4,
                                                                   const Real stiffness) {
    auto *c = new GenericIsometricBendingConstraint();
    const bool res = c->initConstraint(*this, particle1, particle2, particle3, particle4, stiffness);
    if (res) m_constraints.push_back(c);
    return res;
}

bool GenericConstraintsModel::addGenericHingeJoint(const unsigned int rbIndex1,
                                                   const unsigned int rbIndex2,
                                                   const Vector3r &pos,
                                                   const Vector3r &axis) {
    auto *c = new GenericHingeJoint();
    const bool res = c->initConstraint(*this, rbIndex1, rbIndex2, pos, axis);
    if (res) m_constraints.push_back(c);
    return res;
}

bool GenericConstraintsModel::addGenericSliderJoint(const unsigned int rbIndex1,
                                                    const unsigned int rbIndex2,
                                                    const Vector3r &pos,
                                                    const Vector3r &axis) {
    auto *c = new GenericSliderJoint();
    const bool res = c->initConstraint(*this, rbIndex1, rbIndex2, pos, axis);
    if (res) m_constraints.push_back(c);
    return res;
}

bool GenericConstraintsModel::addGenericBallJoint(const unsigned int rbIndex1,
                                                  const unsigned int rbIndex2,
                                                  const Vector3r &pos) {
    auto *c = new GenericBallJoint();
    const bool res = c->initConstraint(*this, rbIndex1, rbIndex2, pos);
    if (res) m_constraints.push_back(c);
    return res;
}