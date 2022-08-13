//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.pbd/line_model.h"

using namespace vox;

line_model::line_model() {
    m_restitutionCoeff = static_cast<Real>(0.6);
    m_frictionCoeff = static_cast<Real>(0.2);
}

line_model::~line_model(void) = default;

line_model::Edges &line_model::getEdges() { return m_edges; }

void line_model::initMesh(const unsigned int nPoints,
                          const unsigned int nQuaternions,
                          const unsigned int indexOffset,
                          const unsigned int indexOffsetQuaternions,
                          unsigned int *indices,
                          unsigned int *indicesQuaternions) {
    m_nPoints = nPoints;
    m_nQuaternions = nQuaternions;
    m_indexOffset = indexOffset;
    m_indexOffsetQuaternions = indexOffsetQuaternions;

    m_edges.resize(nPoints - 1);

    for (unsigned int i = 0; i < nPoints - 1; i++) {
        m_edges[i] = OrientedEdge(indices[2 * i], indices[2 * i + 1], indicesQuaternions[i]);
    }
}

unsigned int line_model::getIndexOffset() const { return m_indexOffset; }

unsigned line_model::getIndexOffsetQuaternions() const { return m_indexOffsetQuaternions; }