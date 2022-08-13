//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <mesh/entity_containers.h>
#include <mesh/triangle_mesh.h>

namespace vox {

FaceIterator FaceContainer::end() const {
  return FaceIterator(static_cast<unsigned int>(m_mesh->nFaces()), m_mesh);
}

FaceConstIterator FaceConstContainer::end() const {
  return FaceConstIterator(static_cast<unsigned int>(m_mesh->nFaces()), m_mesh);
}

VertexIterator VertexContainer::end() const {
  return VertexIterator(static_cast<unsigned int>(m_mesh->nVertices()), m_mesh);
}

VertexConstIterator VertexConstContainer::end() const {
  return VertexConstIterator(static_cast<unsigned int>(m_mesh->nVertices()),
                             m_mesh);
}

} // namespace vox
