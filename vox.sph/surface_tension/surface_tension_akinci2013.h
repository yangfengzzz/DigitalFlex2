//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/surface_tension/surface_tension_base.h"

namespace vox {
/** \brief This class implements the surface tension method introduced
 * by Akinci et al. [ATT13].
 *
 * References:
 * - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids.
 * ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
 */
class SurfaceTension_Akinci2013 : public SurfaceTensionBase {
protected:
    std::vector<Vector3r> m_normals;

public:
    explicit SurfaceTension_Akinci2013(FluidModel *model);
    ~SurfaceTension_Akinci2013() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new SurfaceTension_Akinci2013(model); }

    void step() override;
    void reset() override;

    void computeNormals();

    void performNeighborhoodSearchSort() override;

    FORCE_INLINE Vector3r &getNormal(const unsigned int i) { return m_normals[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getNormal(const unsigned int i) const { return m_normals[i]; }

    FORCE_INLINE void setNormal(const unsigned int i, const Vector3r &val) { m_normals[i] = val; }
};
}  // namespace vox
