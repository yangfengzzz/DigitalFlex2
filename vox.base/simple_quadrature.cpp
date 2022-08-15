//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.base/simple_quadrature.h"

#include "vox.base/logging.h"

using namespace vox;

std::vector<Eigen::Vector3d> SimpleQuadrature::m_samplePoints;
double SimpleQuadrature::m_volume = 0.0;

double SimpleQuadrature::integrate(const Integrand& integrand) {
    double res = 0.0;
    for (auto& m_samplePoint : m_samplePoints) {
        res += m_volume * integrand(m_samplePoint.cast<double>());
    }
    return res;
}

void SimpleQuadrature::determineSamplePointsInSphere(const double radius, unsigned int p) {
    if (p < 1) p = 1;

    m_samplePoints.clear();
    m_samplePoints.reserve(p * p * p);
    const double radius2 = radius * radius;
    const double stepSize = 2.0 * radius / (double)p;
    const double start = -radius + 0.5 * stepSize;
    m_volume = stepSize * stepSize * stepSize;

    Eigen::Vector3d pos;
    pos[0] = start;
    for (unsigned int i = 0; i < p; i++) {
        pos[1] = start;
        for (unsigned int j = 0; j < p; j++) {
            pos[2] = start;
            for (unsigned int k = 0; k < p; k++) {
                // test if sample point is in support radius and if it is not the origin
                const double pn = pos.squaredNorm();
                if (pn < radius2) {
                    m_samplePoints.push_back(pos);
                }
                pos[2] += stepSize;
            }
            pos[1] += stepSize;
        }
        pos[0] += stepSize;
    }
}

void SimpleQuadrature::determineSamplePointsInCircle(const double radius, unsigned int p) {
    if (p < 1) p = 1;

    m_samplePoints.clear();
    m_samplePoints.reserve(p * p);
    const double radius2 = radius * radius;
    const double stepSize = 2.0 * radius / (double)p;
    const double start = -radius + 0.5 * stepSize;
    m_volume = stepSize * stepSize;

    Eigen::Vector3d pos;
    pos[0] = start;
    for (unsigned int i = 0; i < p; i++) {
        pos[1] = start;
        for (unsigned int j = 0; j < p; j++) {
            pos[2] = 0.0;

            // test if sample point is in support radius and if it is not the origin
            const double pn = pos.squaredNorm();
            if (pn < radius2) {
                m_samplePoints.push_back(pos);
            }
            pos[1] += stepSize;
        }
        pos[0] += stepSize;
    }
    LOGI("Number of sampling points: {}", m_samplePoints.size())
}
