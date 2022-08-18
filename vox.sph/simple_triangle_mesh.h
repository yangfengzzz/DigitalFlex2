//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"

namespace vox {
/** \brief Data structure for a triangle mesh with normals and vertex normals.
 */
class SimpleTriangleMesh {
public:
    typedef std::vector<unsigned int> Faces;
    typedef std::vector<Vector3r> Normals;
    typedef std::vector<Vector3r> Vertices;

protected:
    Vertices m_x0;
    Vertices m_x;
    Faces m_indices;
    Normals m_normals;
    Normals m_vertexNormals;

public:
    SimpleTriangleMesh();
    ~SimpleTriangleMesh();

    void release();
    void initMesh(unsigned int nPoints, unsigned int nFaces);
    /** Add a new face.	*/
    void addFace(const unsigned int* indices);
    /** Add a new face.	*/
    void addFace(const int* indices);
    /** Add new vertex. */
    void addVertex(const Vector3r& vertex);

    [[nodiscard]] const Faces& getFaces() const { return m_indices; }
    Faces& getFaces() { return m_indices; }
    [[nodiscard]] const Normals& getFaceNormals() const { return m_normals; }
    Normals& getFaceNormals() { return m_normals; }
    [[nodiscard]] const Normals& getVertexNormals() const { return m_vertexNormals; }
    Normals& getVertexNormals() { return m_vertexNormals; }
    [[nodiscard]] const Vertices& getVertices() const { return m_x; }
    Vertices& getVertices() { return m_x; }
    [[nodiscard]] const Vertices& getVertices0() const { return m_x0; }
    Vertices& getVertices0() { return m_x0; }

    [[nodiscard]] unsigned int numVertices() const { return static_cast<unsigned int>(m_x.size()); }
    [[nodiscard]] unsigned int numFaces() const { return (unsigned int)m_indices.size() / 3; }

    void updateMeshTransformation(const Vector3r& x, const Matrix3r& R);
    void updateNormals();
    void updateVertexNormals();
};

}  // namespace vox