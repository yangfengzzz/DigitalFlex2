//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"
#include "vox.pbd/hashmap.h"

typedef Eigen::Vector3i NeighborhoodSearchCellPos;

namespace vox::utility {
template <>
inline unsigned int hashFunction<NeighborhoodSearchCellPos *>(NeighborhoodSearchCellPos *const &key) {
    const int p1 = 73856093 * (*key)[0];
    const int p2 = 19349663 * (*key)[1];
    const int p3 = 83492791 * (*key)[2];
    return p1 + p2 + p3;
}
}  // namespace vox::utility

namespace vox {
class neighborhood_search_spatial_hashing {
public:
    explicit neighborhood_search_spatial_hashing(unsigned int numParticles = 0,
                                                 Real radius = 0.1,
                                                 unsigned int maxNeighbors = 60u,
                                                 unsigned int maxParticlesPerCell = 50u);
    ~neighborhood_search_spatial_hashing();

    // Spatial hashing
    struct HashEntry {
        HashEntry() = default;
        ;
        unsigned long timestamp{};
        std::vector<unsigned int> particleIndices;
    };

    FORCE_INLINE static int floor(const Real v) {
        return (int)(v + 32768.f) - 32768;  // Shift to get positive values
    }

    void cleanup();
    void neighborhoodSearch(Vector3r *x);
    void neighborhoodSearch(Vector3r *x, unsigned int numBoundaryParticles, Vector3r *boundaryX);
    void update();
    [[nodiscard]] unsigned int **getNeighbors() const;
    [[nodiscard]] unsigned int *getNumNeighbors() const;
    [[nodiscard]] unsigned int getMaxNeighbors() const { return m_maxNeighbors; }

    [[nodiscard]] unsigned int getNumParticles() const;
    void setRadius(Real radius);
    [[nodiscard]] Real getRadius() const;

    [[nodiscard]] FORCE_INLINE unsigned int n_neighbors(unsigned int i) const { return m_numNeighbors[i]; }
    [[nodiscard]] FORCE_INLINE unsigned int neighbor(unsigned int i, unsigned int k) const { return m_neighbors[i][k]; }

private:
    unsigned int m_numParticles;
    unsigned int m_maxNeighbors;
    unsigned int m_maxParticlesPerCell;
    unsigned int **m_neighbors;
    unsigned int *m_numNeighbors;
    Real m_cellGridSize;
    Real m_radius2;
    unsigned int m_currentTimestamp;
    utility::Hashmap<NeighborhoodSearchCellPos *, HashEntry *> m_gridMap;
};
}  // namespace vox