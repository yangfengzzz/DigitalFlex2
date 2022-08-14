//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <string>

#include "vox.base/cubic_lagrange_discrete_grid.h"
#include "vox.pbd/simulation_model.h"
#include "vox.pbd/cubic_sdf_collision_detection.h"
#include "vox.sph/common.h"

namespace vox {
class TimeStepController;

class PBDWrapper {
protected:
    SimulationModel m_model;
    CubicSDFCollisionDetection m_cd;
    TimeStepController *m_timeStep;

    short m_clothSimulationMethod = 2;
    short m_solidSimulationMethod = 2;
    short m_bendingMethod = 2;
    std::string m_sceneName;
    std::string m_sceneFileName;
    bool m_enableMayaExport = false;
    Real m_dampingCoeff = 0.0;

public:
    struct RBData {
        Vector3r x;
        Matrix3r R;
        Vector3r scale;
        std::string objFile;
        int collisionType;
        Real restitution;
        Real friction;
    };

    PBDWrapper();
    ~PBDWrapper();

    void reset();
    void initModel(Real timeStepSize);

    /** Read rigid body scene and create the rigid body model
     */
    void readScene(const std::string &sceneFileName, const std::vector<RBData> &additionalRigidBodies);
    void initTriangleModelConstraints();
    void initTetModelConstraints();

    void timeStep();
    void updateVisModels();

    void loadObj(const std::string &filename, VertexData &vd, utility::IndexedFaceMesh &mesh, const Vector3r &scale);

    SimulationModel &getSimulationModel() { return m_model; }
    DistanceFieldCollisionDetection &getCollisionDetection() { return m_cd; }
    TimeStepController &getTimeStepController();

    [[nodiscard]] Real getDampingCoeff() const { return m_dampingCoeff; }
    void setDampingCoeff(Real val) { m_dampingCoeff = val; }
    [[nodiscard]] short getClothSimulationMethod() const { return m_clothSimulationMethod; }
    void setClothSimulationMethod(short val) { m_clothSimulationMethod = val; }
    [[nodiscard]] short getSolidSimulationMethod() const { return m_solidSimulationMethod; }
    void setSolidSimulationMethod(short val) { m_solidSimulationMethod = val; }
    [[nodiscard]] short getBendingMethod() const { return m_bendingMethod; }
    void setBendingMethod(short val) { m_bendingMethod = val; }
};
}  // namespace vox