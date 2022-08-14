//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"

// ------------------------------------------------------------------------------------
namespace vox {
class MathFunctions {
public:
    /** Implementation of the paper: \n
     * Matthias Müller, Jan Bender, Nuttapong Chentanez and Miles Macklin,
     * "A Robust Method to Extract the Rotational Part of Deformations",
     * ACM SIGGRAPH Motion in Games, 2016
     */
    static void extractRotation(const Matrix3r &A, Quaternionr &q, unsigned int maxIter);

    static void pseudoInverse(const Matrix3r &a, Matrix3r &res);
    static void svdWithInversionHandling(const Matrix3r &A, Vector3r &sigma, Matrix3r &U, Matrix3r &VT);
    static void eigenDecomposition(const Matrix3r &A, Matrix3r &eigenVecs, Vector3r &eigenVals);
    static void jacobiRotate(Matrix3r &A, Matrix3r &R, int p, int q);

    /** Returns two orthogonal vectors to vec which are also orthogonal to each other.
     */
    static void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y);

    /** computes the APD of 8 deformation gradients. (Alg. 3 from the paper: Kugelstadt et al. "Fast Corotated FEM using
     * Operator Splitting", CGF 2018)
     */
    static void APD_Newton(const Matrix3r &F, Quaternionr &q);
};
}  // namespace vox