//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <mesh/entity_iterators.h>
#include <mesh/triangle_mesh.h>

namespace vox {

unsigned int FaceIterator::vertex(unsigned int i) const {
  return m_mesh->faceVertex(m_index, i);
}

FaceIterator::reference FaceIterator::operator*() {
  return m_mesh->face(m_index);
}

FaceConstIterator::reference FaceConstIterator::operator*() {
  return m_mesh->face(m_index);
}

unsigned int &FaceIterator::vertex(unsigned int i) {
  return m_mesh->faceVertex(m_index, i);
}

VertexIterator::reference VertexIterator::operator*() {
  return m_mesh->vertex(m_index);
}

VertexConstIterator::reference VertexConstIterator::operator*() {
  return m_mesh->vertex(m_index);
}

unsigned int VertexIterator::index() const { return m_index; }

IncidentFaceIterator::Self &IncidentFaceIterator::operator++() {
  Halfedge o = m_mesh->opposite(m_h);
  if (o.isBoundary()) {
    m_h = Halfedge();
    return *this;
  }
  m_h = o.next();
  if (m_h == m_begin) {
    m_h = Halfedge();
  }
  return *this;
}

IncidentFaceIterator::IncidentFaceIterator(unsigned int v,
                                           TriangleMesh const *mesh)
    : m_mesh(mesh), m_h(mesh->incident_halfedge(v)),
      m_begin(mesh->incident_halfedge(v)) {
  if (m_h.isBoundary())
    m_h = mesh->opposite(m_h).next();
}

} // namespace vox
