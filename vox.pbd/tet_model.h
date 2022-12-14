//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.base/common.h"
#include "vox.pbd/indexed_face_mesh.h"
#include "vox.pbd/indexed_tet_mesh.h"
#include "vox.pbd/particle_data.h"

namespace vox {
class TetModel {
public:
    TetModel();
    virtual ~TetModel();

    typedef utility::IndexedFaceMesh SurfaceMesh;
    typedef utility::IndexedTetMesh ParticleMesh;

    struct Attachment {
        unsigned int m_index;
        unsigned int m_triIndex;
        Real m_bary[3];
        Real m_dist;
        Real m_minError;
    };

    Vector3r &getInitialX() { return m_initialX; }
    void setInitialX(const Vector3r &val) { m_initialX = val; }
    Matrix3r &getInitialR() { return m_initialR; }
    void setInitialR(const Matrix3r &val) { m_initialR = val; }
    Vector3r &getInitialScale() { return m_initialScale; }
    void setInitialScale(const Vector3r &val) { m_initialScale = val; }

protected:
    /** offset which must be added to get the correct index in the particles array
     */
    unsigned int m_indexOffset{};
    /** Tet mesh of particles which represents the simulation model */
    ParticleMesh m_particleMesh;
    SurfaceMesh m_surfaceMesh;
    VertexData m_visVertices;
    SurfaceMesh m_visMesh;
    Real m_restitutionCoeff;
    Real m_frictionCoeff;
    std::vector<Attachment> m_attachments;
    Vector3r m_initialX;
    Matrix3r m_initialR;
    Vector3r m_initialScale;

    void createSurfaceMesh();
    static void solveQuadraticForZero(const Vector3r &F,
                                      const Vector3r &Fu,
                                      const Vector3r &Fv,
                                      const Vector3r &Fuu,
                                      const Vector3r &Fuv,
                                      const Vector3r &Fvv,
                                      Real &u,
                                      Real &v);
    static bool pointInTriangle(const Vector3r &p0,
                                const Vector3r &p1,
                                const Vector3r &p2,
                                const Vector3r &p,
                                Vector3r &inter,
                                Vector3r &bary);

public:
    SurfaceMesh &getSurfaceMesh();
    VertexData &getVisVertices();
    SurfaceMesh &getVisMesh();
    ParticleMesh &getParticleMesh() { return m_particleMesh; }
    [[nodiscard]] const ParticleMesh &getParticleMesh() const { return m_particleMesh; }
    void cleanupModel();

    [[nodiscard]] unsigned int getIndexOffset() const;

    void initMesh(unsigned int nPoints, unsigned int nTets, unsigned int indexOffset, unsigned int *indices);
    void updateMeshNormals(const ParticleData &pd);

    /** Attach a visualization mesh to the surface of the body.
     * Important: The vertex normals have to be updated before
     * calling this function by calling updateMeshNormals().
     */
    void attachVisMesh(const ParticleData &pd);

    /** Update the visualization mesh of the body.
     * Important: The vertex normals have to be updated before
     * calling this function by calling updateMeshNormals().
     */
    void updateVisMesh(const ParticleData &pd);

    [[nodiscard]] FORCE_INLINE Real getRestitutionCoeff() const { return m_restitutionCoeff; }

    FORCE_INLINE void setRestitutionCoeff(Real val) { m_restitutionCoeff = val; }

    [[nodiscard]] FORCE_INLINE Real getFrictionCoeff() const { return m_frictionCoeff; }

    FORCE_INLINE void setFrictionCoeff(Real val) { m_frictionCoeff = val; }
};
}  // namespace vox