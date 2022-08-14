//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/simulator_base.h"

namespace vox {
/** \brief Base class for data exporters.
 */
class ExporterBase {
protected:
    SimulatorBase* m_base;
    bool m_active;

public:
    explicit ExporterBase(SimulatorBase* base) : m_active(false) { m_base = base; };
    ExporterBase(const ExporterBase&) = delete;
    ExporterBase& operator=(const ExporterBase&) = delete;
    virtual ~ExporterBase() = default;

    virtual void step(unsigned int frame) = 0;
    virtual void reset(){};

    virtual void init(const std::string& outputPath){};

    virtual void setActive(const bool active) { m_active = active; }
    [[nodiscard]] virtual bool getActive() const { return m_active; }
};
}  // namespace vox