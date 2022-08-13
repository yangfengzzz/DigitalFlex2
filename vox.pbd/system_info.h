//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#if WIN32
#define NOMINMAX
#include "windows.h"
#else
#include <limits.h>
#include <unistd.h>
#endif
#include <string>

namespace vox::utility {
class SystemInfo {
public:
    static std::string getHostName() {
#ifdef WIN32
        const unsigned int bufferSize = 32767;
        TCHAR infoBuf[bufferSize];
        DWORD bufCharCount = bufferSize;
        if (!GetComputerName(infoBuf, &bufCharCount)) return "";
        return infoBuf;
#else
        const unsigned int bufferSize = 32767;
        char hostname[bufferSize];
        gethostname(hostname, bufferSize);
        return hostname;
#endif
    }
};
}  // namespace vox::utility