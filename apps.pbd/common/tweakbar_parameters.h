//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <AntTweakBar.h>

#include <vector>

#include "vox.base/common.h"
#include "vox.base/reflect/parameter_object.h"

namespace vox {
class TweakBarParameters {
public:
    typedef std::pair<ParameterObject *, unsigned int> ParameterIndex;

    static void createParameterGUI();
    static void createParameterObjectGUI(ParameterObject *paramObj);

    static void TW_CALL setParameterValue(const void *value, void *clientData);
    static void TW_CALL getParameterValue(void *value, void *clientData);

    static void TW_CALL setTimeStepSizeCB(const void *value, void *clientData);
    static void TW_CALL getTimeStepSizeCB(void *value, void *clientData);

    static void cleanup();

protected:
    static std::vector<std::unique_ptr<ParameterIndex>> m_params;
    static std::vector<std::string> m_objectNames;
};
}  // namespace vox