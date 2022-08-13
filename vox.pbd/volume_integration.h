//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "vox.pbd/common.h"

namespace vox::utility {
class VolumeIntegration {
private:
    int A{};
    int B{};
    int C{};

    // projection integrals
    Real P1{}, Pa{}, Pb{}, Paa{}, Pab{}, Pbb{}, Paaa{}, Paab{}, Pabb{}, Pbbb{};
    // face integrals
    Real Fa{}, Fb{}, Fc{}, Faa{}, Fbb{}, Fcc{}, Faaa{}, Fbbb{}, Fccc{}, Faab{}, Fbbc{}, Fcca{};
    // volume integrals
    Real T0{};
    Real T1[3]{};
    Real T2[3]{};
    Real TP[3]{};

public:
    VolumeIntegration(unsigned int nVertices, unsigned int nFaces, Vector3r *vertices, const unsigned int *indices);

    /** Compute inertia tensor for given geometry and given density.
     */
    void compute_inertia_tensor(Real density);

    /** Return mass of body. */
    [[nodiscard]] Real getMass() const { return m_mass; }
    /** Return volume of body. */
    [[nodiscard]] Real getVolume() const { return m_volume; }
    /** Return inertia tensor of body. */
    [[nodiscard]] Matrix3r const &getInertia() const { return m_theta; }
    /** Return center of mass. */
    [[nodiscard]] Vector3r const &getCenterOfMass() const { return m_r; }

private:
    void volume_integrals();
    void face_integrals(unsigned int i);

    /** Compute various integrations over projection of face.
     */
    void projection_integrals(unsigned int i);

    std::vector<Vector3r> m_face_normals;
    std::vector<Real> m_weights;
    unsigned int m_nVertices;
    unsigned int m_nFaces;
    std::vector<Vector3r> m_vertices;
    const unsigned int *m_indices;

    Real m_mass{}, m_volume{};
    Vector3r m_r;
    Vector3r m_x;
    Matrix3r m_theta;
};
}  // namespace vox::utility