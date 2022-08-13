//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.pbd/neighborhood_search_spatial_hashing.h"

using namespace vox;
using namespace vox::utility;

neighborhood_search_spatial_hashing::neighborhood_search_spatial_hashing(const unsigned int numParticles,
                                                                         const Real radius,
                                                                         const unsigned int maxNeighbors,
                                                                         const unsigned int maxParticlesPerCell)
    : m_gridMap(numParticles * 2) {
    m_cellGridSize = radius;
    m_radius2 = radius * radius;
    m_numParticles = numParticles;
    m_maxParticlesPerCell = maxParticlesPerCell;
    m_maxNeighbors = maxNeighbors;

    m_numNeighbors = nullptr;
    m_neighbors = nullptr;

    if (numParticles != 0) {
        m_numNeighbors = new unsigned int[m_numParticles];
        m_neighbors = new unsigned int *[m_numParticles];
        for (unsigned int i = 0; i < m_numParticles; i++) m_neighbors[i] = new unsigned int[m_maxNeighbors];
    }

    m_currentTimestamp = 0;
}

neighborhood_search_spatial_hashing::~neighborhood_search_spatial_hashing() { cleanup(); }

void neighborhood_search_spatial_hashing::cleanup() {
    for (unsigned int i = 0; i < m_numParticles; i++) delete[] m_neighbors[i];
    delete[] m_neighbors;
    delete[] m_numNeighbors;
    m_numParticles = 0;

    for (unsigned int i = 0; i < m_gridMap.bucket_count(); i++) {
        Hashmap<NeighborhoodSearchCellPos *, neighborhood_search_spatial_hashing::HashEntry *>::KeyValueMap *kvMap =
                m_gridMap.getKeyValueMap(i);
        if (kvMap) {
            for (auto &iter : *kvMap) {
                neighborhood_search_spatial_hashing::HashEntry *entry = iter.second;
                delete entry;
                iter.second = NULL;
            }
        }
    }
}

unsigned int **neighborhood_search_spatial_hashing::getNeighbors() const { return m_neighbors; }

unsigned int *neighborhood_search_spatial_hashing::getNumNeighbors() const { return m_numNeighbors; }

unsigned int neighborhood_search_spatial_hashing::getNumParticles() const { return m_numParticles; }

void neighborhood_search_spatial_hashing::setRadius(const Real radius) {
    m_cellGridSize = radius;
    m_radius2 = radius * radius;
}

Real neighborhood_search_spatial_hashing::getRadius() const { return sqrt(m_radius2); }

void neighborhood_search_spatial_hashing::update() { m_currentTimestamp++; }

void neighborhood_search_spatial_hashing::neighborhoodSearch(Vector3r *x) {
    const Real factor = static_cast<Real>(1.0) / m_cellGridSize;
    for (int i = 0; i < (int)m_numParticles; i++) {
        const int cellPos1 = neighborhood_search_spatial_hashing::floor(x[i][0] * factor) + 1;
        const int cellPos2 = neighborhood_search_spatial_hashing::floor(x[i][1] * factor) + 1;
        const int cellPos3 = neighborhood_search_spatial_hashing::floor(x[i][2] * factor) + 1;
        NeighborhoodSearchCellPos cellPos(cellPos1, cellPos2, cellPos3);
        HashEntry *&entry = m_gridMap[&cellPos];

        if (entry != nullptr) {
            if (entry->timestamp != m_currentTimestamp) {
                entry->timestamp = m_currentTimestamp;
                entry->particleIndices.clear();
            }
        } else {
            auto *newEntry = new HashEntry();
            newEntry->particleIndices.reserve(m_maxParticlesPerCell);
            newEntry->timestamp = m_currentTimestamp;
            entry = newEntry;
        }
        entry->particleIndices.push_back(i);
    }

// loop over all 27 neighboring cells
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)m_numParticles; i++) {
            m_numNeighbors[i] = 0;
            const int cellPos1 = neighborhood_search_spatial_hashing::floor(x[i][0] * factor);
            const int cellPos2 = neighborhood_search_spatial_hashing::floor(x[i][1] * factor);
            const int cellPos3 = neighborhood_search_spatial_hashing::floor(x[i][2] * factor);
            for (unsigned char j = 0; j < 3; j++) {
                for (unsigned char k = 0; k < 3; k++) {
                    for (unsigned char l = 0; l < 3; l++) {
                        NeighborhoodSearchCellPos cellPos(cellPos1 + j, cellPos2 + k, cellPos3 + l);
                        HashEntry *const *entry = m_gridMap.query(&cellPos);

                        if ((entry != nullptr) && (*entry != nullptr) && ((*entry)->timestamp == m_currentTimestamp)) {
                            for (unsigned int pi : (*entry)->particleIndices) {
                                if (pi != i) {
                                    const Real dist2 = (x[i] - x[pi]).squaredNorm();
                                    if (dist2 < m_radius2) {
                                        if (m_numNeighbors[i] < m_maxNeighbors)
                                            m_neighbors[i][m_numNeighbors[i]++] = pi;
                                        //									else
                                        // 											std::cout
                                        // << "too many neighbors detected\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void neighborhood_search_spatial_hashing::neighborhoodSearch(Vector3r *x,
                                                             const unsigned int numBoundaryParticles,
                                                             Vector3r *boundaryX) {
    const Real factor = static_cast<Real>(1.0) / m_cellGridSize;
    for (int i = 0; i < (int)m_numParticles; i++) {
        const int cellPos1 = neighborhood_search_spatial_hashing::floor(x[i][0] * factor) + 1;
        const int cellPos2 = neighborhood_search_spatial_hashing::floor(x[i][1] * factor) + 1;
        const int cellPos3 = neighborhood_search_spatial_hashing::floor(x[i][2] * factor) + 1;
        NeighborhoodSearchCellPos cellPos(cellPos1, cellPos2, cellPos3);
        HashEntry *&entry = m_gridMap[&cellPos];

        if (entry != nullptr) {
            if (entry->timestamp != m_currentTimestamp) {
                entry->timestamp = m_currentTimestamp;
                entry->particleIndices.clear();
            }
        } else {
            auto *newEntry = new HashEntry();
            newEntry->particleIndices.reserve(m_maxParticlesPerCell);
            newEntry->timestamp = m_currentTimestamp;
            entry = newEntry;
        }
        entry->particleIndices.push_back(i);
    }

    for (int i = 0; i < (int)numBoundaryParticles; i++) {
        const int cellPos1 = neighborhood_search_spatial_hashing::floor(boundaryX[i][0] * factor) + 1;
        const int cellPos2 = neighborhood_search_spatial_hashing::floor(boundaryX[i][1] * factor) + 1;
        const int cellPos3 = neighborhood_search_spatial_hashing::floor(boundaryX[i][2] * factor) + 1;
        NeighborhoodSearchCellPos cellPos(cellPos1, cellPos2, cellPos3);
        HashEntry *&entry = m_gridMap[&cellPos];

        if (entry != nullptr) {
            if (entry->timestamp != m_currentTimestamp) {
                entry->timestamp = m_currentTimestamp;
                entry->particleIndices.clear();
            }
        } else {
            auto *newEntry = new HashEntry();
            newEntry->particleIndices.reserve(m_maxParticlesPerCell);
            newEntry->timestamp = m_currentTimestamp;
            entry = newEntry;
        }
        entry->particleIndices.push_back(m_numParticles + i);
    }

// loop over all 27 neighboring cells
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)m_numParticles; i++) {
            m_numNeighbors[i] = 0;
            const int cellPos1 = neighborhood_search_spatial_hashing::floor(x[i][0] * factor);
            const int cellPos2 = neighborhood_search_spatial_hashing::floor(x[i][1] * factor);
            const int cellPos3 = neighborhood_search_spatial_hashing::floor(x[i][2] * factor);
            for (unsigned char j = 0; j < 3; j++) {
                for (unsigned char k = 0; k < 3; k++) {
                    for (unsigned char l = 0; l < 3; l++) {
                        NeighborhoodSearchCellPos cellPos(cellPos1 + j, cellPos2 + k, cellPos3 + l);
                        HashEntry *const *entry = m_gridMap.query(&cellPos);

                        if ((entry != nullptr) && (*entry != nullptr) && ((*entry)->timestamp == m_currentTimestamp)) {
                            for (unsigned int pi : (*entry)->particleIndices) {
                                if (pi != i) {
                                    Real dist2;
                                    if (pi < m_numParticles)
                                        dist2 = (x[i] - x[pi]).squaredNorm();
                                    else
                                        dist2 = (x[i] - boundaryX[pi - m_numParticles]).squaredNorm();

                                    if (dist2 < m_radius2) {
                                        if (m_numNeighbors[i] < m_maxNeighbors)
                                            m_neighbors[i][m_numNeighbors[i]++] = pi;
                                        // else
                                        // std::cout << "too many neighbors detected\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}