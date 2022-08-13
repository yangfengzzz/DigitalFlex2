//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <array>
#include <iterator>
#include <vector>

#include "vox.pbd/common.h"

namespace vox::utility {
class IndexedFaceMesh {
public:
    struct Edge {
        std::array<unsigned int, 2> m_face;
        std::array<unsigned int, 2> m_vert;
    };

public:
    typedef std::vector<unsigned int> Faces;
    typedef std::vector<Vector3r> FaceNormals;
    typedef std::vector<Vector3r> VertexNormals;
    typedef std::vector<std::vector<unsigned int>> FacesEdges;
    typedef std::vector<Edge> Edges;
    typedef std::vector<std::vector<unsigned int>> VerticesEdges;
    typedef std::vector<std::vector<unsigned int>> VerticesFaces;
    typedef std::vector<unsigned int> UVIndices;
    typedef std::vector<Vector2r> UVs;

protected:
    unsigned int m_numPoints{};
    Faces m_indices;
    Edges m_edges;
    FacesEdges m_facesEdges;
    bool m_closed{};
    UVIndices m_uvIndices;
    UVs m_uvs;
    VerticesFaces m_verticesFaces;
    VerticesEdges m_verticesEdges;
    const unsigned int m_verticesPerFace = 3u;
    FaceNormals m_normals;
    VertexNormals m_vertexNormals;
    bool m_flatShading{};

public:
    IndexedFaceMesh();
    IndexedFaceMesh(IndexedFaceMesh const &other);
    IndexedFaceMesh &operator=(IndexedFaceMesh const &other);
    ~IndexedFaceMesh();

    void release();
    [[nodiscard]] bool isClosed() const;
    [[nodiscard]] bool getFlatShading() const { return m_flatShading; }
    void setFlatShading(const bool v) { m_flatShading = v; }
    void initMesh(unsigned int nPoints, unsigned int nEdges, unsigned int nFaces);
    void addFace(const unsigned int *indices);
    void addFace(const int *indices);
    void addUV(Real u, Real v);
    void addUVIndex(unsigned int index);

    [[nodiscard]] const Faces &getFaces() const { return m_indices; }
    Faces &getFaces() { return m_indices; }
    [[nodiscard]] const FaceNormals &getFaceNormals() const { return m_normals; }
    FaceNormals &getFaceNormals() { return m_normals; }
    [[nodiscard]] const VertexNormals &getVertexNormals() const { return m_vertexNormals; }
    VertexNormals &getVertexNormals() { return m_vertexNormals; }
    Edges &getEdges() { return m_edges; }
    [[nodiscard]] const Edges &getEdges() const { return m_edges; }
    [[nodiscard]] const FacesEdges &getFacesEdges() const { return m_facesEdges; }
    [[nodiscard]] const UVIndices &getUVIndices() const { return m_uvIndices; }
    [[nodiscard]] const UVs &getUVs() const { return m_uvs; }
    [[nodiscard]] const VerticesFaces &getVertexFaces() const { return m_verticesFaces; }
    [[nodiscard]] const VerticesEdges &getVertexEdges() const { return m_verticesEdges; }

    [[nodiscard]] unsigned int numVertices() const { return m_numPoints; }
    [[nodiscard]] unsigned int numFaces() const { return (unsigned int)m_indices.size() / m_verticesPerFace; }
    [[nodiscard]] unsigned int numEdges() const { return (unsigned int)m_edges.size(); }
    [[nodiscard]] unsigned int numUVs() const { return (unsigned int)m_uvs.size(); }

    void copyUVs(const UVIndices &uvIndices, const UVs &uvs);

    void buildNeighbors();

    template <class PositionData>
    void updateNormals(const PositionData &pd, unsigned int offset);

    template <class PositionData>
    void updateVertexNormals(const PositionData &pd);

    [[nodiscard]] unsigned int getVerticesPerFace() const;
};

template <class PositionData>
void IndexedFaceMesh::updateNormals(const PositionData &pd, const unsigned int offset) {
    m_normals.resize(numFaces());

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numFaces(); i++) {
            // Get first three points of face
            const Vector3r &a = pd.getPosition(m_indices[m_verticesPerFace * i] + offset);
            const Vector3r &b = pd.getPosition(m_indices[m_verticesPerFace * i + 1] + offset);
            const Vector3r &c = pd.getPosition(m_indices[m_verticesPerFace * i + 2] + offset);

            // Create normal
            Vector3r v1 = b - a;
            Vector3r v2 = c - a;

            m_normals[i] = v1.cross(v2);
            m_normals[i].normalize();
            // fix normals of degenerate triangles that can become zero vectors
            if (m_normals[i].squaredNorm() < 1e-6f) m_normals[i] = Vector3r::UnitX();
        }
    }
}

template <class PositionData>
void IndexedFaceMesh::updateVertexNormals(const PositionData &pd) {
    m_vertexNormals.resize(numVertices());

    for (unsigned int i = 0; i < numVertices(); i++) {
        m_vertexNormals[i].setZero();
    }

    for (unsigned int i = 0u; i < numFaces(); i++) {
        const Vector3r &n = m_normals[i];
        m_vertexNormals[m_indices[m_verticesPerFace * i]] += n;
        m_vertexNormals[m_indices[m_verticesPerFace * i + 1]] += n;
        m_vertexNormals[m_indices[m_verticesPerFace * i + 2]] += n;
    }

    for (unsigned int i = 0; i < numVertices(); i++) {
        m_vertexNormals[i].normalize();
    }
}

}  // namespace vox::utility
