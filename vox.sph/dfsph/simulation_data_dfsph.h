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
/** \brief Simulation data which is required by the method Divergence-free Smoothed Particle Hydrodynamics introduced
 * by Bender and Koschier [BK15,BK17].
 *
 * References:
 * - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics
 * Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL:
 * http://doi.acm.org/10.1145/2786784.2786796
 * - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on
 * Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
 */
class SimulationDataDFSPH {
public:
    SimulationDataDFSPH();
    virtual ~SimulationDataDFSPH();

protected:
    /** \brief factor \f$\alpha_i\f$ */
    std::vector<std::vector<Real>> m_factor;
    /** \brief stores \f$\kappa\f$ value of last time step for a warm start of the pressure solver */
    std::vector<std::vector<Real>> m_kappa;
    /** \brief stores \f$\kappa^v\f$ value of last time step for a warm start of the divergence solver */
    std::vector<std::vector<Real>> m_kappaV;
    /** \brief advected density */
    std::vector<std::vector<Real>> m_density_adv;

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
    void emittedParticles(FluidModel* model, unsigned int startIndex);

    [[nodiscard]] FORCE_INLINE Real getFactor(const unsigned int fluidIndex, const unsigned int i) const {
        return m_factor[fluidIndex][i];
    }

    FORCE_INLINE Real& getFactor(const unsigned int fluidIndex, const unsigned int i) {
        return m_factor[fluidIndex][i];
    }

    FORCE_INLINE void setFactor(const unsigned int fluidIndex, const unsigned int i, const Real p) {
        m_factor[fluidIndex][i] = p;
    }

    [[nodiscard]] FORCE_INLINE Real getKappa(const unsigned int fluidIndex, const unsigned int i) const {
        return m_kappa[fluidIndex][i];
    }

    FORCE_INLINE Real& getKappa(const unsigned int fluidIndex, const unsigned int i) { return m_kappa[fluidIndex][i]; }

    FORCE_INLINE void setKappa(const unsigned int fluidIndex, const unsigned int i, const Real p) {
        m_kappa[fluidIndex][i] = p;
    }

    [[nodiscard]] FORCE_INLINE Real getKappaV(const unsigned int fluidIndex, const unsigned int i) const {
        return m_kappaV[fluidIndex][i];
    }

    FORCE_INLINE Real& getKappaV(const unsigned int fluidIndex, const unsigned int i) {
        return m_kappaV[fluidIndex][i];
    }

    FORCE_INLINE void setKappaV(const unsigned int fluidIndex, const unsigned int i, const Real p) {
        m_kappaV[fluidIndex][i] = p;
    }

    [[nodiscard]] FORCE_INLINE Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const {
        return m_density_adv[fluidIndex][i];
    }

    FORCE_INLINE Real& getDensityAdv(const unsigned int fluidIndex, const unsigned int i) {
        return m_density_adv[fluidIndex][i];
    }

    FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d) {
        m_density_adv[fluidIndex][i] = d;
    }
};
}  // namespace vox