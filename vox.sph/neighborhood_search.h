//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

// #define GPU_NEIGHBORHOOD_SEARCH
#if USE_DOUBLE
#define CUNSEARCH_USE_DOUBLE_PRECISION
#endif

#ifdef GPU_NEIGHBORHOOD_SEARCH
#include "cuNSearch.h"
typedef cuNSearch::NeighborhoodSearch NeighborhoodSearch;
#else
#include "vox.base/compact_search.h"
typedef CompactNSearch::NeighborhoodSearch NeighborhoodSearch;
#endif
