//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/mesh/entity_iterators.h"

namespace vox {

class TriangleMesh;

class FaceContainer {
public:
    [[nodiscard]] FaceIterator begin() const { return FaceIterator(0, m_mesh); }
    [[nodiscard]] FaceIterator end() const;

private:
    friend class TriangleMesh;
    FaceContainer(TriangleMesh *mesh) : m_mesh(mesh) {}

    TriangleMesh *m_mesh;
};

class FaceConstContainer {
public:
    [[nodiscard]] FaceConstIterator begin() const { return FaceConstIterator(0, m_mesh); }
    [[nodiscard]] FaceConstIterator end() const;

private:
    friend class TriangleMesh;
    FaceConstContainer(TriangleMesh const *mesh) : m_mesh(mesh) {}

    TriangleMesh const *m_mesh;
};

class IncidentFaceContainer {
public:
    [[nodiscard]] IncidentFaceIterator begin() const { return IncidentFaceIterator(m_v, m_mesh); }
    [[nodiscard]] IncidentFaceIterator end() const { return IncidentFaceIterator(); }

private:
    friend class TriangleMesh;
    IncidentFaceContainer(unsigned int v, TriangleMesh const *mesh) : m_v(v), m_mesh(mesh) {}

    TriangleMesh const *m_mesh;
    unsigned int m_v;
};

class VertexContainer {
public:
    [[nodiscard]] VertexIterator begin() const { return {0, m_mesh}; }
    [[nodiscard]] VertexIterator end() const;

private:
    friend class TriangleMesh;
    VertexContainer(TriangleMesh *mesh) : m_mesh(mesh) {}

    TriangleMesh *m_mesh;
};

class VertexConstContainer {
public:
    [[nodiscard]] VertexConstIterator begin() const { return VertexConstIterator(0, m_mesh); }
    [[nodiscard]] VertexConstIterator end() const;

private:
    friend class TriangleMesh;
    VertexConstContainer(TriangleMesh const *mesh) : m_mesh(mesh) {}

    TriangleMesh const *m_mesh;
};
}  // namespace vox
