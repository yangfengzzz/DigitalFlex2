//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "apps.pbd/common/mini_gl.h"
#include "vox.base/common.h"

namespace vox {
class Visualization {
public:
    template <class PositionData>
    static void drawMesh(const PositionData &pd,
                         const utility::IndexedFaceMesh &mesh,
                         unsigned int offset,
                         const float *color);
    template <class PositionData>
    static void drawTexturedMesh(const PositionData &pd,
                                 const utility::IndexedFaceMesh &mesh,
                                 unsigned int offset,
                                 const float *color);
};

template <class PositionData>
void Visualization::drawMesh(const PositionData &pd,
                             const utility::IndexedFaceMesh &mesh,
                             const unsigned int offset,
                             const float *const color) {
    // draw mesh
    const Vector3r *vertexNormals = mesh.getVertexNormals().data();

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &pd.getPosition(offset)[0]);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_REAL, GL_FALSE, 0, &vertexNormals[0][0]);
    } else {
        float speccolor[4] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f);
        glColor3fv(color);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_REAL, 0, &pd.getPosition(0)[0]);
        glNormalPointer(GL_REAL, 0, &vertexNormals[0][0]);
    }

    glDrawElements(GL_TRIANGLES, (GLsizei)3 * mesh.numFaces(), GL_UNSIGNED_INT, mesh.getFaces().data());

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(2);
    } else {
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
}

template <class PositionData>
void Visualization::drawTexturedMesh(const PositionData &pd,
                                     const utility::IndexedFaceMesh &mesh,
                                     const unsigned int offset,
                                     const float *const color) {
    // draw mesh
    const Vector3r *vertexNormals = mesh.getVertexNormals().data();
    const Vector2r *uvs = mesh.getUVs().data();

    MiniGL::bindTexture();

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &pd.getPosition(offset)[0]);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_REAL, GL_FALSE, 0, &uvs[0][0]);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_REAL, GL_FALSE, 0, &vertexNormals[0][0]);
    } else {
        float speccolor[4] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
        glColor3fv(color);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glVertexPointer(3, GL_REAL, 0, &pd.getPosition(0)[0]);
        glTexCoordPointer(2, GL_REAL, 0, &uvs[0][0]);
        glNormalPointer(GL_REAL, 0, &vertexNormals[0][0]);
    }

    glDrawElements(GL_TRIANGLES, (GLsizei)3 * mesh.numFaces(), GL_UNSIGNED_INT, mesh.getFaces().data());

    if (MiniGL::checkOpenGLVersion(3, 3)) {
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
    } else {
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    }

    MiniGL::unbindTexture();
}
}  // namespace vox