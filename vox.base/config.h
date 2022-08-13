//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

namespace CompactNSearch {
#ifdef USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif
}  // namespace CompactNSearch

#define INITIAL_NUMBER_OF_INDICES 50
#define INITIAL_NUMBER_OF_NEIGHBORS 50
