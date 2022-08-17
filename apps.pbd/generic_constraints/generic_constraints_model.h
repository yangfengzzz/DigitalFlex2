//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.pbd/indexed_face_mesh.h"
#include "vox.pbd/particle_data.h"
#include "vox.pbd/simulation_model.h"

namespace vox {
class GenericConstraintsModel : public SimulationModel {
public:
    GenericConstraintsModel();
    ~GenericConstraintsModel() override;

    bool addGenericDistanceConstraint(unsigned int particle1, unsigned int particle2, Real stiffness);
    bool addGenericIsometricBendingConstraint(unsigned int particle1,
                                              unsigned int particle2,
                                              unsigned int particle3,
                                              unsigned int particle4,
                                              Real stiffness);
    bool addGenericHingeJoint(unsigned int rbIndex1, unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis);
    bool addGenericSliderJoint(unsigned int rbIndex1, unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis);
    bool addGenericBallJoint(unsigned int rbIndex1, unsigned int rbIndex2, const Vector3r &pos);
};
}  // namespace vox