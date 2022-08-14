//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <random>
#include <string>
#include <unordered_map>

#include "vox.base/common.h"

namespace vox {
/** \brief This class implements a Poisson disk sampling for the surface
 * of 3D models.
 */
class PoissonDiskSampling {
    typedef Eigen::Matrix<int, 3, 1, Eigen::DontAlign> CellPos;

    struct CellPosHasher {
        std::size_t operator()(const CellPos &k) const {
            const int p1 = 73856093 * k[0];
            const int p2 = 19349663 * k[1];
            const int p3 = 83492791 * k[2];
            return (size_t)(p1 + p2 + p3);
        }
    };

public:
    PoissonDiskSampling();

    /** \brief Struct to store the information of the initial points
     */
    struct InitialPointInfo {
        CellPos cP;
        Vector3r pos;
        unsigned int ID;
    };

    /** \brief Struct to store the hash entry (spatial hashing)
     */
    struct HashEntry {
        HashEntry() { startIndex = 0; };
        std::vector<unsigned int> samples;
        unsigned int startIndex;
    };

    FORCE_INLINE static int floor(const Real v) {
        return (int)(v + 32768.f) - 32768;  // Shift to get positive values
    }

    /** Performs the poisson sampling with the
     * respective parameters. Compare
     * http://graphics.cs.umass.edu/pubs/sa_2010.pdf
     *
     * @param mesh mesh data of sampled body
     * @param vertices vertex data of sampled data
     * @param sampledVertices sampled vertices that will be returned
     * @param minRadius minimal distance of sampled vertices
     * @param numTestpointsPerFace # of generated test points per face of body
     * @param distanceNorm 0: euclidean norm, 1: approx geodesic distance
     * @param numTrials # of iterations used to find samples
     */
    void sampleMesh(unsigned int numVertices,
                    const Vector3r *vertices,
                    unsigned int numFaces,
                    const unsigned int *faces,
                    Real minRadius,
                    unsigned int numTrials,
                    unsigned int distanceNorm,
                    std::vector<Vector3r> &samples);

private:
    Real m_r{};
    unsigned int m_numTrials{};
    unsigned int m_numTestpointsPerFace{};
    unsigned int m_distanceNorm{};
    std::vector<Vector3r> m_faceNormals;
    std::vector<Real> m_areas;
    Real m_totalArea{};

    Real m_cellSize{};
    Vector3r m_minVec;
    Vector3r m_maxVec;

    std::vector<InitialPointInfo> m_initialInfoVec;
    std::vector<std::vector<CellPos>> m_phaseGroups;

    Real m_maxArea{};

    void computeFaceNormals(unsigned int numVertices,
                            const Vector3r *vertices,
                            unsigned int numFaces,
                            const unsigned int *faces);
    void determineTriangleAreas(unsigned int numVertices,
                                const Vector3r *vertices,
                                unsigned int numFaces,
                                const unsigned int *faces);
    void generateInitialPointSet(unsigned int numVertices,
                                 const Vector3r *vertices,
                                 unsigned int numFaces,
                                 const unsigned int *faces);
    unsigned int getAreaIndex(const std::vector<Real> &areas,
                              Real totalArea,
                              std::default_random_engine &generator,
                              std::uniform_real_distribution<Real> &distribution) const;
    void parallelUniformSurfaceSampling(std::vector<Vector3r> &samples);

    void quickSort(int left, int right);
    int partition(int left, int right);
    static bool compareCellID(CellPos &a, CellPos &b);

    void determineMinX(unsigned int numVertices, const Vector3r *vertices);

    bool nbhConflict(const std::unordered_map<CellPos, HashEntry, CellPosHasher> &kvMap, const InitialPointInfo &iPI);
    bool checkCell(const std::unordered_map<CellPos, HashEntry, CellPosHasher> &kvMap,
                   const CellPos &cell,
                   const InitialPointInfo &iPI);
};
}  // namespace vox