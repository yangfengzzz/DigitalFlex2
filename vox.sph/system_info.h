//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/common.h"
#if WIN32
#define NOMINMAX
#include "windows.h"
#else
#include <unistd.h>

#include <climits>
#endif

namespace vox::utility {
class SystemInfo {
public:
    static std::string getHostName() {
#ifdef WIN32
        const unsigned int bufferSize = 32767;
        TCHAR *infoBuf = new TCHAR[bufferSize];
        DWORD bufCharCount = bufferSize;
        if (!GetComputerName(infoBuf, &bufCharCount)) return "";
        std::string name = infoBuf;
        delete[] infoBuf;
        return name;
#else
        const unsigned int bufferSize = 32767;
        char hostname[bufferSize];
        gethostname(hostname, bufferSize);
        return hostname;
#endif
    }
};
}  // namespace vox::utility