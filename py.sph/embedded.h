//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <pybind11/pybind11.h>

#include <functional>
#include <iostream>
#include <string>

namespace vox {
class SimulatorBase;

class Embedded {
private:
    static Embedded* current;

protected:
    pybind11::object m_mainScope;
    std::vector<std::string> m_functions;

    void findFunctions(const std::string& moduleName);

public:
    Embedded();
    ~Embedded();

    // Singleton
    static Embedded* getCurrent();
    static void setCurrent(Embedded* tm);
    static bool hasCurrent();

    std::string import_script(const std::string& fileName);
    void exec_script_str(const std::string& str);
    void exec_script_file(const std::string& fileName);
    void exec_fct(const std::string& moduleName, const std::string& fctName);
    void exec_init(const std::string& moduleName, SimulatorBase* base);

    [[nodiscard]] const std::vector<std::string>& getFunctions() const { return m_functions; }

    void reloadModule(const std::string& moduleName);
    void releaseModule(const std::string& moduleName);
};
}  // namespace vox