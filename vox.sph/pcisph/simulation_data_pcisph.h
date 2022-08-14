//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Simulation data which is required by the method Predictive-corrective Incompressible SPH introduced
 * by Solenthaler and Pajarola [SP09].
 *
 * References:
 * - [SP09] B. Solenthaler and R. Pajarola. Predictive-corrective incompressible SPH. ACM Trans. Graph.,
 * 28(3):40:1-40:6, July 2009. URL: http://doi.acm.org/10.1145/1531326.1531346
 */
class SimulationDataPCISPH {
public:
    SimulationDataPCISPH();
    virtual ~SimulationDataPCISPH();

protected:
    std::vector<Real> m_pcisph_factor;

    std::vector<std::vector<Vector3r>> m_predX;
    std::vector<std::vector<Vector3r>> m_predV;
    std::vector<std::vector<Real>> m_densityAdv;
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

    Real getPCISPH_ScalingFactor(const unsigned int fluidIndex) { return m_pcisph_factor[fluidIndex]; }

    void emittedParticles(FluidModel *model, unsigned int startIndex);

    FORCE_INLINE Vector3r &getPredictedPosition(const unsigned int fluidIndex, const unsigned int i) {
        return m_predX[fluidIndex][i];
    }

    [[nodiscard]] FORCE_INLINE const Vector3r &getPredictedPosition(const unsigned int fluidIndex,
                                                                    const unsigned int i) const {
        return m_predX[fluidIndex][i];
    }

    FORCE_INLINE void setPredictedPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r &pos) {
        m_predX[fluidIndex][i] = pos;
    }

    FORCE_INLINE Vector3r &getPredictedVelocity(const unsigned int fluidIndex, const unsigned int i) {
        return m_predV[fluidIndex][i];
    }

    [[nodiscard]] FORCE_INLINE const Vector3r &getPredictedVelocity(const unsigned int fluidIndex,
                                                                    const unsigned int i) const {
        return m_predV[fluidIndex][i];
    }

    FORCE_INLINE void setPredictedVelocity(const unsigned int fluidIndex, const unsigned int i, const Vector3r &vel) {
        m_predV[fluidIndex][i] = vel;
    }

    [[nodiscard]] FORCE_INLINE Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const {
        return m_densityAdv[fluidIndex][i];
    }

    FORCE_INLINE Real &getDensityAdv(const unsigned int fluidIndex, const unsigned int i) {
        return m_densityAdv[fluidIndex][i];
    }

    FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d) {
        m_densityAdv[fluidIndex][i] = d;
    }

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