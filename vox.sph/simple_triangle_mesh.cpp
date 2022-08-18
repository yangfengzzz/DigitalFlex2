//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/simple_triangle_mesh.h"

using namespace vox;

SimpleTriangleMesh::SimpleTriangleMesh() = default;

SimpleTriangleMesh::~SimpleTriangleMesh() { release(); }

void SimpleTriangleMesh::initMesh(const unsigned int nPoints, const unsigned int nFaces) {
    m_x0.reserve(nPoints);
    m_x.reserve(nPoints);
    m_indices.reserve(nFaces * 3);
    m_normals.reserve(nFaces);
    m_vertexNormals.reserve(nPoints);
}

void SimpleTriangleMesh::release() {
    m_indices.clear();
    m_x0.clear();
    m_x.clear();
    m_normals.clear();
    m_vertexNormals.clear();
}

void SimpleTriangleMesh::addFace(const unsigned int *const indices) {
    for (unsigned int i = 0u; i < 3; i++) m_indices.push_back(indices[i]);
}

void SimpleTriangleMesh::addFace(const int *const indices) {
    for (unsigned int i = 0u; i < 3; i++) m_indices.push_back((unsigned int)indices[i]);
}

void SimpleTriangleMesh::addVertex(const Vector3r &vertex) {
    m_x0.push_back(vertex);
    m_x.push_back(vertex);
}

void SimpleTriangleMesh::updateNormals() {
    m_normals.resize(numFaces());

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numFaces(); i++) {
            // Get first three points of face
            const Vector3r &a = m_x[m_indices[3 * i]];
            const Vector3r &b = m_x[m_indices[3 * i + 1]];
            const Vector3r &c = m_x[m_indices[3 * i + 2]];

            // Create normal
            Vector3r v1 = b - a;
            Vector3r v2 = c - a;

            m_normals[i] = v1.cross(v2);
            m_normals[i].normalize();
        }
    }
}

void SimpleTriangleMesh::updateVertexNormals() {
    m_vertexNormals.resize(numVertices());

    for (unsigned int i = 0; i < numVertices(); i++) {
        m_vertexNormals[i].setZero();
    }

    for (unsigned int i = 0u; i < numFaces(); i++) {
        const Vector3r &n = m_normals[i];
        m_vertexNormals[m_indices[3 * i]] += n;
        m_vertexNormals[m_indices[3 * i + 1]] += n;
        m_vertexNormals[m_indices[3 * i + 2]] += n;
    }

    for (unsigned int i = 0; i < numVertices(); i++) {
        m_vertexNormals[i].normalize();
    }
}

void SimpleTriangleMesh::updateMeshTransformation(const Vector3r &x, const Matrix3r &R) {
    for (unsigned int j = 0; j < numVertices(); j++) m_x[j] = R * m_x0[j] + x;
}
