//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <array>
#include <vector>

#include "vox.base/common.h"

namespace vox::utility {
class IndexedTetMesh {
public:
    struct Edge {
        std::array<unsigned int, 2> m_vert;
    };

    struct Face {
        // edge indices
        std::array<unsigned int, 3> m_edges;
        // tet indices
        std::array<unsigned int, 2> m_tets;
    };

    struct Tet {
        std::array<unsigned int, 6> m_edges;
        std::array<unsigned int, 4> m_faces;
    };

public:
    typedef std::vector<unsigned int> Tets;
    typedef std::vector<unsigned int> Faces;
    typedef std::vector<Tet> TetData;
    typedef std::vector<Face> FaceData;
    typedef std::vector<Edge> Edges;
    typedef std::vector<std::vector<unsigned int>> VerticesTets;
    typedef std::vector<std::vector<unsigned int>> VerticesFaces;
    typedef std::vector<std::vector<unsigned int>> VerticesEdges;

protected:
    unsigned int m_numPoints{};
    Tets m_tetIndices;
    Faces m_faceIndices;
    Edges m_edges;
    FaceData m_faces;
    TetData m_tets;
    VerticesTets m_verticesTets;
    VerticesFaces m_verticesFaces;
    VerticesEdges m_verticesEdges;

public:
    IndexedTetMesh();
    ~IndexedTetMesh();

    void release();
    void initMesh(unsigned int nPoints, unsigned int nEdges, unsigned int nFaces, unsigned int nTets);
    void addTet(const unsigned int *indices);
    void addTet(const int *indices);

    [[nodiscard]] const Faces &getFaces() const { return m_faceIndices; }
    Faces &getFaces() { return m_faceIndices; }
    [[nodiscard]] const Tets &getTets() const { return m_tetIndices; }
    Tets &getTets() { return m_tetIndices; }
    Edges &getEdges() { return m_edges; }
    [[nodiscard]] const Edges &getEdges() const { return m_edges; }
    [[nodiscard]] const FaceData &getFaceData() const { return m_faces; }
    [[nodiscard]] const TetData &getTetData() const { return m_tets; }
    [[nodiscard]] const VerticesTets &getVertexTets() const { return m_verticesTets; }
    [[nodiscard]] const VerticesFaces &getVertexFaces() const { return m_verticesFaces; }
    [[nodiscard]] const VerticesEdges &getVertexEdges() const { return m_verticesEdges; }

    [[nodiscard]] unsigned int numVertices() const { return m_numPoints; }
    [[nodiscard]] unsigned int numFaces() const { return (unsigned int)m_faceIndices.size() / 3; }
    [[nodiscard]] unsigned int numTets() const { return (unsigned int)m_tetIndices.size() / 4; }
    [[nodiscard]] unsigned int numEdges() const { return (unsigned int)m_edges.size(); }

    void buildNeighbors();
};
}  // namespace vox::utility