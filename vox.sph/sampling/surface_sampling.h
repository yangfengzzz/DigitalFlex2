//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/utilities/poisson_disk_sampling.h"
#include "vox.sph/utilities/regular_sampling_2d.h"
#include "vox.sph/utilities/regular_triangle_sampling.h"

namespace vox {
enum SurfaceSamplingMode { PoissonDisk, RegularTriangle, Regular2D };
}