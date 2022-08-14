//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <string>

#include "vox.base/reflect/parameter_object.h"

namespace vox {
class SimulatorBase;

class ScriptObject : public ParameterObject {
protected:
    SimulatorBase* m_base;
    bool m_scriptLoaded;
    std::string m_scriptModule;
    std::string m_scriptFile;

    void initParameters() override;
    void addFunctionParameters();
    void removeFunctionParameters();

public:
    static int SCRIPT_FILE;

    explicit ScriptObject(SimulatorBase* base);
    ~ScriptObject() override;

    void init();
    std::string loadScriptFile(const std::string& fileName);

    void execResetFct();
    void execStepFct();

    void updateFunctionParameters();
};
}  // namespace vox