//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/discrete_grid/discrete_grid.h"
#include "vox.sph/boundary_model.h"
#include "vox.base/common.h"
#include "vox.sph/sph_kernels.h"

namespace vox {
class TimeStep;

/** \brief The boundary model stores the information required for boundary handling
 * using the approach of Bender et al. 2019 [BKWK19].
 *
 * References:
 * - [BKWK19] Jan Bender, Tassilo Kugelstadt, Marcel Weiler, and Dan Koschier. Volume maps: an implicit boundary
 * representation for SPH. In Proceedings of ACM SIGGRAPH Conference on Motion, Interaction and Games, MIG '19. ACM,
 * 2019. URL: https://dl.acm.org/doi/10.1145/3359566.3360077
 */
class BoundaryModel_Bender2019 : public BoundaryModel {
public:
    BoundaryModel_Bender2019();
    ~BoundaryModel_Bender2019() override;

protected:
    // Density or volume map
    DiscreteGrid *m_map;
    // values required for volume maps
    std::vector<std::vector<Real>> m_boundaryVolume;
    std::vector<std::vector<Vector3r>> m_boundaryXj;
    // maxmimal distance of a mesh point to the center of mass (required for CFL)
    Real m_maxDist;
    Real m_maxVel;

public:
    void initModel(RigidBodyObject *rbo);

    void reset() override;

    DiscreteGrid *getMap() { return m_map; }
    void setMap(DiscreteGrid *map) { m_map = map; }

    [[nodiscard]] Real getMaxDist() const { return m_maxDist; }
    void setMaxDist(Real val) { m_maxDist = val; }

    [[nodiscard]] Real getMaxVel() const { return m_maxVel; }
    void setMaxVel(Real val) { m_maxVel = val; }

    [[nodiscard]] FORCE_INLINE const Real &getBoundaryVolume(const unsigned int fluidIndex,
                                                             const unsigned int i) const {
        return m_boundaryVolume[fluidIndex][i];
    }

    FORCE_INLINE Real &getBoundaryVolume(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryVolume[fluidIndex][i];
    }

    FORCE_INLINE void setBoundaryVolume(const unsigned int fluidIndex, const unsigned int i, const Real &val) {
        m_boundaryVolume[fluidIndex][i] = val;
    }

    FORCE_INLINE Vector3r &getBoundaryXj(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryXj[fluidIndex][i];
    }

    [[nodiscard]] FORCE_INLINE const Vector3r &getBoundaryXj(const unsigned int fluidIndex,
                                                             const unsigned int i) const {
        return m_boundaryXj[fluidIndex][i];
    }

    FORCE_INLINE void setBoundaryXj(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val) {
        m_boundaryXj[fluidIndex][i] = val;
    }
};
}  // namespace vox