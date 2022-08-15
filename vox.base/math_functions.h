//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"

// ------------------------------------------------------------------------------------
namespace vox {
class MathFunctions {
public:
    static void jacobiRotate(Matrix3r &A, Matrix3r &R, int p, int q);

    /**
     * Return the inf norm of the matrix.
     */
    static Real infNorm(const Matrix3r &A);
    /**
     * Return the one norm of the matrix.
     */
    static Real oneNorm(const Matrix3r &A);
    
    /**
     * Implementation of the paper: \n
     * Matthias M�ller, Jan Bender, Nuttapong Chentanez and Miles Macklin,
     * "A Robust Method to Extract the Rotational Part of Deformations",
     * ACM SIGGRAPH Motion in Games, 2016
     */
    static void extractRotation(const Matrix3r &A, Quaternionr &q, unsigned int maxIter);
    
    static void pseudoInverse(const Matrix3r &a, Matrix3r &res);

    /**
     * Perform a singular value decomposition of matrix A: A = U * sigma * V^T.
     * This function returns two proper rotation matrices U and V^T which do not
     * contain a reflection. Reflections are corrected by the inversion handling
     * proposed by Irving et al. 2004.
     */
    static void svdWithInversionHandling(const Matrix3r &A, Vector3r &sigma, Matrix3r &U, Matrix3r &VT);
    
    static void eigenDecomposition(const Matrix3r &A, Matrix3r &eigenVecs, Vector3r &eigenVals);

    /**
     * Returns two orthogonal vectors to vec which are also orthogonal to each other.
     */
    static void getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y);

    /**
     * computes the APD of 8 deformation gradients. (Alg. 3 from the paper: Kugelstadt et al. "Fast Corotated FEM using
     * Operator Splitting", CGF 2018)
     */
    static void APD_Newton(const Matrix3r &F, Quaternionr &q);
    
    /**
     * Perform polar decomposition A = (U D U^T) R
     */
    static void polarDecomposition(const Matrix3r &A, Matrix3r &R, Matrix3r &U, Matrix3r &D);

    /**
     * Perform a polar decomposition of matrix M and return the rotation matrix R.
     * This method handles the degenerated cases.
     */
    static void polarDecompositionStable(const Matrix3r &M, Real tolerance, Matrix3r &R);

    static Real cotTheta(const Vector3r &v, const Vector3r &w);

    /** Computes the cross product matrix of a vector.
     * @param  v		input vector
     * @param  v_hat	resulting cross product matrix
     */
    static void crossProductMatrix(const Vector3r &v, Matrix3r &v_hat);


};
}  // namespace vox
