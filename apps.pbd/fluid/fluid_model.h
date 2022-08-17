//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.pbd/neighborhood_search_spatial_hashing.h"
#include "vox.pbd/particle_data.h"

namespace vox {
class FluidModel {
public:
    FluidModel();
    virtual ~FluidModel();

protected:
    Real viscosity;
    Real m_density0;
    Real m_particleRadius;
    Real m_supportRadius{};
    ParticleData m_particles;
    std::vector<Vector3r> m_boundaryX;
    std::vector<Real> m_boundaryPsi;
    std::vector<Real> m_density;
    std::vector<Real> m_lambda;
    std::vector<Vector3r> m_deltaX;
    NeighborhoodSearchSpatialHashing* m_neighborhoodSearch;

    void initMasses();

    void resizeFluidParticles(unsigned int newSize);
    void releaseFluidParticles();

public:
    void cleanupModel();
    virtual void reset();

    ParticleData& getParticles();

    void initModel(unsigned int nFluidParticles,
                   Vector3r* fluidParticles,
                   unsigned int nBoundaryParticles,
                   Vector3r* boundaryParticles);

    [[nodiscard]] unsigned int numBoundaryParticles() const { return (unsigned int)m_boundaryX.size(); }
    [[nodiscard]] Real getDensity0() const { return m_density0; }
    [[nodiscard]] Real getSupportRadius() const { return m_supportRadius; }
    [[nodiscard]] Real getParticleRadius() const { return m_particleRadius; }
    void setParticleRadius(Real val) {
        m_particleRadius = val;
        m_supportRadius = static_cast<Real>(4.0) * m_particleRadius;
    }
    NeighborhoodSearchSpatialHashing* getNeighborhoodSearch() { return m_neighborhoodSearch; }

    [[nodiscard]] Real getViscosity() const { return viscosity; }
    void setViscosity(Real val) { viscosity = val; }

    [[nodiscard]] FORCE_INLINE const Vector3r& getBoundaryX(const unsigned int i) const { return m_boundaryX[i]; }

    FORCE_INLINE Vector3r& getBoundaryX(const unsigned int i) { return m_boundaryX[i]; }

    FORCE_INLINE void setBoundaryX(const unsigned int i, const Vector3r& val) { m_boundaryX[i] = val; }

    [[nodiscard]] FORCE_INLINE const Real& getBoundaryPsi(const unsigned int i) const { return m_boundaryPsi[i]; }

    FORCE_INLINE Real& getBoundaryPsi(const unsigned int i) { return m_boundaryPsi[i]; }

    FORCE_INLINE void setBoundaryPsi(const unsigned int i, const Real& val) { m_boundaryPsi[i] = val; }

    [[nodiscard]] FORCE_INLINE const Real& getLambda(const unsigned int i) const { return m_lambda[i]; }

    FORCE_INLINE Real& getLambda(const unsigned int i) { return m_lambda[i]; }

    FORCE_INLINE void setLambda(const unsigned int i, const Real& val) { m_lambda[i] = val; }

    [[nodiscard]] FORCE_INLINE const Real& getDensity(const unsigned int i) const { return m_density[i]; }

    FORCE_INLINE Real& getDensity(const unsigned int i) { return m_density[i]; }

    FORCE_INLINE void setDensity(const unsigned int i, const Real& val) { m_density[i] = val; }

    FORCE_INLINE Vector3r& getDeltaX(const unsigned int i) { return m_deltaX[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r& getDeltaX(const unsigned int i) const { return m_deltaX[i]; }

    FORCE_INLINE void setDeltaX(const unsigned int i, const Vector3r& val) { m_deltaX[i] = val; }
};
}  // namespace vox