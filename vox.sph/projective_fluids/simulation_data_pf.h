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
/** \brief Simulation data which is required by the method Projective Fluids introduced by Weiler, Koschier
 * and Bender [WKB16].
 *
 * References:
 * - [WKB16] Marcel Weiler, Dan Koschier, and Jan Bender. Projective fluids. In Proceedings of the 9th International
 * Conference on Motion in Games, MIG '16, 79-84. New York, NY, USA, 2016. ACM. URL:
 * http://doi.acm.org/10.1145/2994258.2994282
 */
class SimulationDataPF {
public:
    SimulationDataPF();
    virtual ~SimulationDataPF();

protected:
    /** \brief particle position from last timestep */
    std::vector<std::vector<Vector3r>> m_old_position;

    /** \brief number of neighbors that are fluid particles */
    std::vector<std::vector<unsigned int>> m_num_fluid_neighbors;

    /** \brief positions predicted from momentum */
    std::vector<std::vector<Vector3r>> m_s;

    /** \brief diagonal of system matrix, used by preconditioner*/
    std::vector<std::vector<Vector3r>> m_mat_diag;

    std::vector<unsigned int> m_particleOffset;

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

    [[nodiscard]] FORCE_INLINE Vector3r getOldPosition(const unsigned int fluidIndex, const unsigned int i) const {
        return m_old_position[fluidIndex][i];
    }

    FORCE_INLINE Vector3r& getOldPosition(const unsigned int fluidIndex, const unsigned int i) {
        return m_old_position[fluidIndex][i];
    }

    FORCE_INLINE void setOldPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r& p) {
        m_old_position[fluidIndex][i] = p;
    }

    [[nodiscard]] FORCE_INLINE unsigned int getNumFluidNeighbors(const unsigned int fluidIndex,
                                                                 const unsigned int i) const {
        return m_num_fluid_neighbors[fluidIndex][i];
    }

    FORCE_INLINE unsigned int& getNumFluidNeighbors(const unsigned int fluidIndex, const unsigned int i) {
        return m_num_fluid_neighbors[fluidIndex][i];
    }

    FORCE_INLINE void setNumFluidNeighbors(const unsigned int fluidIndex, const unsigned int i, const unsigned int n) {
        m_num_fluid_neighbors[fluidIndex][i] = n;
    }

    [[nodiscard]] FORCE_INLINE const Vector3r& getS(const unsigned int fluidIndex, const unsigned int i) const {
        return m_s[fluidIndex][i];
    }

    FORCE_INLINE Vector3r& getS(const unsigned int fluidIndex, const unsigned int i) { return m_s[fluidIndex][i]; }

    FORCE_INLINE void setS(const unsigned int fluidIndex, const unsigned int i, const Vector3r& s) {
        m_s[fluidIndex][i] = s;
    }

    [[nodiscard]] FORCE_INLINE const Vector3r& getDiag(const unsigned int fluidIndex, const unsigned int i) const {
        return m_mat_diag[fluidIndex][i];
    }

    FORCE_INLINE Vector3r& getDiag(const unsigned int fluidIndex, const unsigned int i) {
        return m_mat_diag[fluidIndex][i];
    }

    FORCE_INLINE void setDiag(const unsigned int fluidIndex, const unsigned int i, const Vector3r& s) {
        m_mat_diag[fluidIndex][i] = s;
    }

    [[nodiscard]] FORCE_INLINE const unsigned int& getParticleOffset(const unsigned int fluidIndex) const {
        return m_particleOffset[fluidIndex];
    }
};
}  // namespace vox