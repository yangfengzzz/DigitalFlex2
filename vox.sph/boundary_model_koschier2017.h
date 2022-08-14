//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/discrete_grid/discrete_grid.h"
#include "vox.sph/boundary_model.h"
#include "vox.sph/common.h"
#include "vox.sph/sph_kernels.h"

namespace vox {
class TimeStep;

/** \brief The boundary model stores the information required for boundary handling
 * using the approach of Koschier and Bender 2017 [KB17].
 *
 * References:
 * - [KB17] Dan Koschier and Jan Bender. Density maps for improved SPH boundary handling. In ACM SIGGRAPH/Eurographics
 * Symposium on Computer Animation, 1-10. July 2017. URL: http://dx.doi.org/10.1145/3099564.3099565
 */
class BoundaryModel_Koschier2017 : public BoundaryModel {
public:
    BoundaryModel_Koschier2017();
    ~BoundaryModel_Koschier2017() override;

protected:
    // Density map
    DiscreteGrid *m_map;
    // values required for density maps
    std::vector<std::vector<Real>> m_boundaryDensity;
    std::vector<std::vector<Vector3r>> m_boundaryDensityGradient;
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

    [[nodiscard]] FORCE_INLINE const Real &getBoundaryDensity(const unsigned int fluidIndex,
                                                              const unsigned int i) const {
        return m_boundaryDensity[fluidIndex][i];
    }

    FORCE_INLINE Real &getBoundaryDensity(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryDensity[fluidIndex][i];
    }

    FORCE_INLINE void setBoundaryDensity(const unsigned int fluidIndex, const unsigned int i, const Real &val) {
        m_boundaryDensity[fluidIndex][i] = val;
    }

    FORCE_INLINE Vector3r &getBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryDensityGradient[fluidIndex][i];
    }

    [[nodiscard]] FORCE_INLINE const Vector3r &getBoundaryDensityGradient(const unsigned int fluidIndex,
                                                                          const unsigned int i) const {
        return m_boundaryDensityGradient[fluidIndex][i];
    }

    FORCE_INLINE void setBoundaryDensityGradient(const unsigned int fluidIndex,
                                                 const unsigned int i,
                                                 const Vector3r &val) {
        m_boundaryDensityGradient[fluidIndex][i] = val;
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