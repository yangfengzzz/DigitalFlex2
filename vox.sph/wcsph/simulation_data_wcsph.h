//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Simulation data which is required by the method Weakly Compressible SPH for Free Surface Flows introduced
 * by Becker and Teschner [BT07].
 *
 * References:
 * - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM
 * SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland,
 * 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
 */
class SimulationDataWCSPH {
public:
    SimulationDataWCSPH();
    virtual ~SimulationDataWCSPH();

protected:
    std::vector<std::vector<Real>> m_pressure;
    std::vector<std::vector<Vector3r>> m_pressureAccel;

public:
    /** Initialize the arrays containing the particle data.
     */
    virtual void init();

    /** Release the arrays containing the particle data.
     */
    virtual void cleanup();

    /** Reset the particle data.
     */
    virtual void reset();

    /** Important: First call m_model->performNeighborhoodSearchSort()
     * to call the z_sort of the neighborhood search.
     */
    void performNeighborhoodSearchSort();

    void emittedParticles(FluidModel *model, unsigned int startIndex);

    [[nodiscard]] FORCE_INLINE Real getPressure(const unsigned int fluidIndex, const unsigned int i) const {
        return m_pressure[fluidIndex][i];
    }

    FORCE_INLINE Real &getPressure(const unsigned int fluidIndex, const unsigned int i) {
        return m_pressure[fluidIndex][i];
    }

    FORCE_INLINE void setPressure(const unsigned int fluidIndex, const unsigned int i, const Real p) {
        m_pressure[fluidIndex][i] = p;
    }

    FORCE_INLINE Vector3r &getPressureAccel(const unsigned int fluidIndex, const unsigned int i) {
        return m_pressureAccel[fluidIndex][i];
    }

    [[nodiscard]] FORCE_INLINE const Vector3r &getPressureAccel(const unsigned int fluidIndex,
                                                                const unsigned int i) const {
        return m_pressureAccel[fluidIndex][i];
    }

    FORCE_INLINE void setPressureAccel(const unsigned int fluidIndex, const unsigned int i, const Vector3r &val) {
        m_pressureAccel[fluidIndex][i] = val;
    }
};
}  // namespace vox