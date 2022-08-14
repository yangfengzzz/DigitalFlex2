//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/animation_field_system.h"

#include "vox.sph/simulation.h"

using namespace vox;

AnimationFieldSystem::AnimationFieldSystem() : m_fields() {}

AnimationFieldSystem::~AnimationFieldSystem() {
    for (auto &m_field : m_fields) delete m_field;
}

void AnimationFieldSystem::step() {
    for (auto &m_field : m_fields) {
        m_field->step();
    }
}

void AnimationFieldSystem::reset() {
    for (auto &m_field : m_fields) {
        m_field->reset();
    }
}

void AnimationFieldSystem::addAnimationField(const std::string &particleFieldName,
                                             const Vector3r &pos,
                                             const Matrix3r &rotation,
                                             const Vector3r &scale,
                                             const std::string expression[3],
                                             const unsigned int type) {
    m_fields.push_back(new AnimationField(particleFieldName, pos, rotation, scale, expression, type));
}
