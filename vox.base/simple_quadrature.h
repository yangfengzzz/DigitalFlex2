//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "vox.base/common.h"

namespace vox {
class SimpleQuadrature {
public:
    using Integrand = std::function<double(Eigen::Vector3d const&)>;
    using Domain = Eigen::AlignedBox3d;

    static std::vector<Eigen::Vector3d> m_samplePoints;
    static double m_volume;

    static void determineSamplePointsInSphere(double radius, unsigned int p);
    static void determineSamplePointsInCircle(double radius, unsigned int p);
    static double integrate(const Integrand& integrand);
};
}  // namespace vox