//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <iostream>
#include <unordered_map>

#include "vox.base/logging.h"
#include "vox.sph/common.h"
#include "vox.sph/timing.h"

namespace vox::utility {
#define INCREASE_COUNTER(counterName, increaseBy) vox::utility::Counting::increaseCounter(counterName, increaseBy);

#define INIT_COUNTING \
    std::unordered_map<std::string, vox::utility::AverageCount> vox::utility::Counting::m_averageCounts;

struct AverageCount {
    Real sum;
    unsigned int numberOfCalls;
};

class Counting {
public:
    static std::unordered_map<std::string, AverageCount> m_averageCounts;

    static void reset() { m_averageCounts.clear(); }

    FORCE_INLINE static void increaseCounter(const std::string &name, const Real increaseBy) {
        std::unordered_map<std::string, AverageCount>::iterator iter;
        iter = Counting::m_averageCounts.find(name);
        if (iter != Counting::m_averageCounts.end()) {
            iter->second.sum += increaseBy;
            iter->second.numberOfCalls++;
        } else {
            AverageCount ac{};
            ac.sum = increaseBy;
            ac.numberOfCalls = 1;
            Counting::m_averageCounts[name] = ac;
        }
    }

    FORCE_INLINE static void printAverageCounts() {
        std::unordered_map<std::string, AverageCount>::iterator iter;
        for (iter = Counting::m_averageCounts.begin(); iter != Counting::m_averageCounts.end(); iter++) {
            AverageCount &ac = iter->second;
            const double avgCount = ac.sum / ac.numberOfCalls;
            LOGI("Average number: {}:{}", iter->first.c_str(), avgCount)
        }
        LOGI("-----------------------------------------------------")
    }

    FORCE_INLINE static void printCounterSums() {
        std::unordered_map<std::string, AverageCount>::iterator iter;
        for (iter = Counting::m_averageCounts.begin(); iter != Counting::m_averageCounts.end(); iter++) {
            AverageCount &ac = iter->second;
            LOGI("Total number: {}:{}", iter->first.c_str(), ac.sum)
        }
        LOGI("-----------------------------------------------------")
    }
};
}  // namespace vox::utility