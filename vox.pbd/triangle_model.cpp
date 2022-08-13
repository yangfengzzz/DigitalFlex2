//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "triangle_model.h"

#include "vox.pbd/position_based_dynamics.h"
#include "vox.pbd/position_based_rigid_body_dynamics.h"

using namespace vox;

TriangleModel::TriangleModel() : m_particleMesh() {
    m_restitutionCoeff = static_cast<Real>(0.6);
    m_frictionCoeff = static_cast<Real>(0.2);
}

TriangleModel::~TriangleModel(void) { cleanupModel(); }

void TriangleModel::cleanupModel() { m_particleMesh.release(); }

void TriangleModel::updateMeshNormals(const ParticleData &pd) {
    m_particleMesh.updateNormals(pd, m_indexOffset);
    m_particleMesh.updateVertexNormals(pd);
}

void TriangleModel::initMesh(const unsigned int nPoints,
                             const unsigned int nFaces,
                             const unsigned int indexOffset,
                             unsigned int *indices,
                             const ParticleMesh::UVIndices &uvIndices,
                             const ParticleMesh::UVs &uvs) {
    m_indexOffset = indexOffset;
    m_particleMesh.release();

    m_particleMesh.initMesh(nPoints, nFaces * 2, nFaces);

    for (unsigned int i = 0; i < nFaces; i++) {
        m_particleMesh.addFace(&indices[3 * i]);
    }
    m_particleMesh.copyUVs(uvIndices, uvs);
    m_particleMesh.buildNeighbors();
}

unsigned int TriangleModel::getIndexOffset() const { return m_indexOffset; }
