//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.pbd/common.h"
#include "vox.pbd/indexed_face_mesh.h"
#include "vox.pbd/particle_data.h"

namespace vox {
class RigidBodyGeometry {
public:
    RigidBodyGeometry();
    virtual ~RigidBodyGeometry();

    typedef Utilities::IndexedFaceMesh Mesh;

protected:
    Mesh m_mesh;
    VertexData m_vertexData_local;
    VertexData m_vertexData;

public:
    Mesh &getMesh();
    VertexData &getVertexData();
    const VertexData &getVertexData() const;
    VertexData &getVertexDataLocal();
    const VertexData &getVertexDataLocal() const;

    void initMesh(const unsigned int nVertices,
                  const unsigned int nFaces,
                  const Vector3r *vertices,
                  const unsigned int *indices,
                  const Mesh::UVIndices &uvIndices,
                  const Mesh::UVs &uvs,
                  const Vector3r &scale = Vector3r(1.0, 1.0, 1.0),
                  const bool flatShading = false);
    void updateMeshTransformation(const Vector3r &x, const Matrix3r &R);
    void updateMeshNormals(const VertexData &vd);
};
}  // namespace vox