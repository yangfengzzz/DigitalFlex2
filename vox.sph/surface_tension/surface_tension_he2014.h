//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/surface_tension/surface_tension_base.h"

namespace vox {
/** \brief This class implements the surface tension method introduced
 * by He et al. [HWZ+14].
 *
 * References:
 * - [HWZ+14] Xiaowei He, Huamin Wang, Fengjun Zhang, Hongan Wang, Guoping Wang, and Kun Zhou. Robust simulation of
 * sparsely sampled thin features in SPH-based free surface flows. ACM Trans. Graph., 34(1):7:1-7:9, December 2014. URL:
 * http://doi.acm.org/10.1145/2682630
 */
class SurfaceTension_He2014 : public SurfaceTensionBase {
protected:
    std::vector<Real> m_color;
    std::vector<Real> m_gradC2;

public:
    explicit SurfaceTension_He2014(FluidModel* model);
    ~SurfaceTension_He2014() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_He2014(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    [[nodiscard]] FORCE_INLINE Real getColor(const unsigned int i) const { return m_color[i]; }

    FORCE_INLINE Real& getColor(const unsigned int i) { return m_color[i]; }

    FORCE_INLINE void setColor(const unsigned int i, const Real p) { m_color[i] = p; }

    [[nodiscard]] FORCE_INLINE Real getGradC2(const unsigned int i) const { return m_gradC2[i]; }

    FORCE_INLINE Real& getGradC2(const unsigned int i) { return m_gradC2[i]; }

    FORCE_INLINE void setGradC2(const unsigned int i, const Real p) { m_gradC2[i] = p; }
};
}  // namespace vox