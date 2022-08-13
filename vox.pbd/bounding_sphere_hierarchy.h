//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.pbd/bounding_sphere.h"
#include "vox.pbd/common.h"
#include "vox.pbd/kdtree.h"

namespace vox {
class PointCloudBSH : public KDTree<bounding_sphere> {
public:
    using super = KDTree<bounding_sphere>;

    PointCloudBSH();

    void init(const Vector3r *vertices, const unsigned int numVertices);
    Vector3r const &entity_position(unsigned int i) const final;
    void compute_hull(unsigned int b, unsigned int n, bounding_sphere &hull) const final;
    void compute_hull_approx(unsigned int b, unsigned int n, bounding_sphere &hull) const final;

private:
    const Vector3r *m_vertices;
    unsigned int m_numVertices;
};

class TetMeshBSH : public KDTree<bounding_sphere> {
public:
    using super = KDTree<bounding_sphere>;

    TetMeshBSH();

    void init(const Vector3r *vertices,
              const unsigned int numVertices,
              const unsigned int *indices,
              const unsigned int numTets,
              const Real tolerance);
    Vector3r const &entity_position(unsigned int i) const final;
    void compute_hull(unsigned int b, unsigned int n, bounding_sphere &hull) const final;
    void compute_hull_approx(unsigned int b, unsigned int n, bounding_sphere &hull) const final;
    void updateVertices(const Vector3r *vertices);

private:
    const Vector3r *m_vertices;
    unsigned int m_numVertices;
    const unsigned int *m_indices;
    unsigned int m_numTets;
    Real m_tolerance;
    std::vector<Vector3r> m_com;
};

class BVHTest {
public:
    using TraversalCallback = std::function<void(unsigned int node_index1, unsigned int node_index2)>;

    static void traverse(PointCloudBSH const &b1, TetMeshBSH const &b2, TraversalCallback func);
    static void traverse(PointCloudBSH const &b1,
                         const unsigned int node_index1,
                         TetMeshBSH const &b2,
                         const unsigned int node_index2,
                         TraversalCallback func);
};
}  // namespace vox