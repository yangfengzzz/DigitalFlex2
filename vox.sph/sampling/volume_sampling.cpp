//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/sampling/volume_sampling.h"

#include "vox.base/sdf_functions.h"
#include "vox.base/timing.h"

using namespace Eigen;
using namespace std;
namespace vox::utility {

void VolumeSampling::sampleMesh(const unsigned int numVertices,
                                const Vector3r *vertices,
                                const unsigned int numFaces,
                                const unsigned int *faces,
                                const Real radius,
                                const AlignedBox3r *region,
                                const std::array<unsigned int, 3> &resolution,
                                const bool invert,
                                const unsigned int sampleMode,
                                std::vector<Vector3r> &samples) {
    AlignedBox3r bbox = SDFFunctions::computeBoundingBox(numVertices, vertices);

    if (region) {
        for (unsigned int i = 0; i < 3; i++) {
            bbox.min()[i] = std::max(region->min()[i], bbox.min()[i]);
            bbox.max()[i] = std::min(region->max()[i], bbox.max()[i]);
        }
    }

    CubicLagrangeDiscreteGrid *sdf =
            SDFFunctions::generateSDF(numVertices, vertices, numFaces, faces, bbox, resolution, invert);

    const Real diameter = static_cast<Real>(2.0) * radius;

    // sample object
    const unsigned int numberOfSamplePoints =
            (((unsigned int)((1.0f / diameter) * (bbox.max()[2] - bbox.min()[2]))) + 1) *
            (((unsigned int)((1.0f / diameter) * (bbox.max()[1] - bbox.min()[1]))) + 1) *
            (((unsigned int)((1.0f / diameter) * (bbox.max()[0] - bbox.min()[0]))) + 1);

    unsigned int currentSample = 0;
    //	Real currentPercent = 0.0;
    int counter_x = 0;
    int counter_y = 0;
    Real xshift = diameter;
    Real yshift = diameter;

    if (sampleMode == 1)
        yshift = sqrt(static_cast<Real>(3.0)) * radius;
    else if (sampleMode == 2) {
        xshift = sqrt(static_cast<Real>(3.0)) * radius;
        yshift = sqrt(static_cast<Real>(6.0)) * diameter / static_cast<Real>(3.0);
    }
    for (Real z = bbox.min()[2]; z <= bbox.max()[2]; z += diameter) {
        for (Real y = bbox.min()[1]; y <= bbox.max()[1]; y += yshift) {
            for (Real x = bbox.min()[0]; x <= bbox.max()[0]; x += xshift) {
                Vector3r particlePosition;
                if (sampleMode == 1) {
                    if (counter_y % 2 == 0)
                        particlePosition = Vector3r(x, y + radius, z + radius);
                    else
                        particlePosition = Vector3r(x + radius, y + radius, z);
                } else if (sampleMode == 2) {
                    particlePosition = Vector3r(x, y + radius, z + radius);

                    Vector3r shift_vec(0, 0, 0);
                    if (counter_x % 2) {
                        shift_vec[2] += diameter / (static_cast<Real>(2.0) * (counter_y % 2 ? -1 : 1));
                    }
                    if (counter_y % 2) {
                        shift_vec[0] += xshift / static_cast<Real>(2.0);
                        shift_vec[2] += diameter / static_cast<Real>(2.0);
                    }
                    particlePosition += shift_vec;
                } else {
                    // Use center of voxel
                    particlePosition = Vector3r(x + radius, y + radius, z + radius);
                }

                if (SDFFunctions::distance(sdf, particlePosition, 0.0) < 0.0) samples.push_back(particlePosition);
                currentSample++;

                //  				if ((Real)currentSample / (Real)numberOfSamplePoints > currentPercent)
                //  				{
                //  					LOG_INFO << currentPercent * 100.0 << "%";
                //  					currentPercent += 0.05;
                //  				}
                counter_x++;
            }
            counter_x = 0;
            counter_y++;
        }
        counter_y = 0;
    }
    // 	LOG_INFO << "100%";

    delete sdf;
}
}  // namespace vox::utility
