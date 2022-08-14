//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Eigen/Dense>

namespace vox {
class GaussQuadrature {
public:
    using Integrand = std::function<double(Eigen::Vector3d const&)>;
    using Domain = Eigen::AlignedBox3d;

    static double integrate(Integrand integrand, Domain const& domain, unsigned int p);
    static void exportSamples(unsigned int p);
};
}  // namespace vox