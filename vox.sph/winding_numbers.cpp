//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/winding_numbers.h"

#include <cmath>

using namespace vox::utility;
using namespace vox;
using namespace Eigen;

Real WindingNumbers::computeGeneralizedWindingNumber(const Vector3r& p_,
                                                     const Vector3r& a_,
                                                     const Vector3r& b_,
                                                     const Vector3r& c_) {
    const Vector3r& a = a_ - p_;
    const Vector3r& b = b_ - p_;
    const Vector3r& c = c_ - p_;

    const Real normA = a.norm();
    const Real normB = b.norm();
    const Real normC = c.norm();

    Matrix3r A;
    A.row(0) = a;
    A.row(1) = b;
    A.row(2) = c;
    const Real det = A.determinant();
    const Real divisor = normA * normB * normC + (a.dot(b)) * normC + (b.dot(c)) * normA + (c.dot(a)) * normB;

    static const Real tau = 2.0 * M_PI;
    return std::atan2(det, divisor) / tau;  // Only divide by 2*pi instead of 4*pi because there was a 2 out front
}

Real WindingNumbers::computeGeneralizedWindingNumber(const Vector3r& p, const SimpleTriangleMesh& mesh) {
    const unsigned int* faces = mesh.getFaces().data();
    const unsigned int nFaces = mesh.numFaces();
    const Vector3r* v = mesh.getVertices().data();

    // compute generalized winding number
    Real w_p = 0;
#pragma omp parallel for reduction(+ : w_p)
    for (int idx = 0; idx < static_cast<int>(nFaces); idx++) {
        const Vector3r& a = v[faces[idx * 3 + 0]];
        const Vector3r& b = v[faces[idx * 3 + 1]];
        const Vector3r& c = v[faces[idx * 3 + 2]];

        w_p += computeGeneralizedWindingNumber(p, a, b, c);
    }

    return w_p;
}
