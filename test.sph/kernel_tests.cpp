//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <gtest/gtest.h>

#include "vox.sph/common.h"
#include "vox.sph/sph_kernels.h"

using namespace vox;

#ifdef USE_DOUBLE
double eps = 1.0e-5;
#else
float eps = 1.0e-4f;
#endif

template <typename TestType>
void test3DKernel() {
    const Real supportRadius = static_cast<Real>(4.0 * 0.025);
    TestType::setRadius(supportRadius);
    const unsigned int numberOfSteps = 50;
    const Real stepSize = static_cast<Real>(2.0) * supportRadius / (Real)(numberOfSteps - 1);
    Vector3r xi;
    xi.setZero();
    Real sum = 0.0;
    Vector3r sumV = Vector3r::Zero();
    bool positive = true;
    Real V = std::pow(stepSize, 3.f);
    for (unsigned int i = 0; i < numberOfSteps; i++) {
        for (unsigned int j = 0; j < numberOfSteps; j++) {
            for (unsigned int k = 0; k < numberOfSteps; k++) {
                const Vector3r xj(-supportRadius + i * stepSize, -supportRadius + j * stepSize,
                                  -supportRadius + k * stepSize);
                const Real W = TestType::W(xi - xj);
                sum += W * V;
                sumV += TestType::gradW(xi - xj) * V;
                if (W < -eps) positive = false;
            }
        }
    }
    ASSERT_TRUE(fabs(sum - 1.0) < eps);
    ASSERT_TRUE(sumV.norm() < eps);
    ASSERT_TRUE(positive);
}

TEST(SPHKERNEL, ThreeD) {
    test3DKernel<CubicKernel>();
    test3DKernel<Poly6Kernel>();
    test3DKernel<SpikyKernel>();
    test3DKernel<WendlandQuinticC2Kernel>();
    test3DKernel<PrecomputedKernel<CubicKernel>>();
}

template <typename TestType>
void test2DKernel() {
    const Real supportRadius = static_cast<Real>(4.0 * 0.025);
    TestType::setRadius(supportRadius);
    const unsigned int numberOfSteps = 50;
    const Real stepSize = static_cast<Real>(2.0) * supportRadius / (Real)(numberOfSteps - 1);
    Vector3r xi;
    xi.setZero();
    Real sum = 0.0;
    Vector3r sumV = Vector3r::Zero();
    bool positive = true;
    Real V = std::pow(stepSize, 2.f);
    for (unsigned int i = 0; i < numberOfSteps; i++) {
        for (unsigned int j = 0; j < numberOfSteps; j++) {
            const Vector3r xj(-supportRadius + i * stepSize, -supportRadius + j * stepSize, 0.0);
            const Real W = TestType::W(xi - xj);
            sum += W * V;
            sumV += TestType::gradW(xi - xj) * V;
            if (W < -eps) positive = false;
        }
    }
    ASSERT_TRUE(fabs(sum - 1.0) < eps);
    ASSERT_TRUE(sumV.norm() < eps);
    ASSERT_TRUE(positive);
}

TEST(SPHKERNEL, TWOD) {
    test2DKernel<CubicKernel2D>();
    test2DKernel<WendlandQuinticC2Kernel2D>();
    test2DKernel<PrecomputedKernel<CubicKernel2D>>();
}