//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <memory>

#include "vox.base/cubic_lagrange_discrete_grid.h"
#include "vox.pbd/common.h"
#include "vox.pbd/distance_field_collision_detection.h"

namespace vox {
/** Collision detection based on cubic signed distance fields.
 */
class CubicSDFCollisionDetection : public DistanceFieldCollisionDetection {
public:
    using Grid = CubicLagrangeDiscreteGrid;
    using GridPtr = std::shared_ptr<CubicLagrangeDiscreteGrid>;

    struct CubicSDFCollisionObject : public DistanceFieldCollisionDetection::DistanceFieldCollisionObject {
        std::string m_sdfFile;
        Vector3r m_scale;
        GridPtr m_sdf;
        static int TYPE_ID;

        CubicSDFCollisionObject();
        ~CubicSDFCollisionObject() override;
        [[nodiscard]] int &getTypeId() const override { return TYPE_ID; }
        bool collisionTest(
                const Vector3r &x, Real tolerance, Vector3r &cp, Vector3r &n, Real &dist, Real maxDist = 0.0) override;
        double distance(const Eigen::Vector3d &x, Real tolerance) override;
    };

public:
    CubicSDFCollisionDetection();
    ~CubicSDFCollisionDetection() override;

    bool isDistanceFieldCollisionObject(CollisionObject *co) const override;

    void addCubicSDFCollisionObject(unsigned int bodyIndex,
                                    unsigned int bodyType,
                                    const Vector3r *vertices,
                                    unsigned int numVertices,
                                    const std::string &sdfFile,
                                    const Vector3r &scale,
                                    bool testMesh = true,
                                    bool invertSDF = false);
    void addCubicSDFCollisionObject(unsigned int bodyIndex,
                                    unsigned int bodyType,
                                    const Vector3r *vertices,
                                    unsigned int numVertices,
                                    GridPtr sdf,
                                    const Vector3r &scale,
                                    bool testMesh = true,
                                    bool invertSDF = false);
};
}  // namespace vox