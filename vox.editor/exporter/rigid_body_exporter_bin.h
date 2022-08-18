//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.editor/exporter/exporter_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Rigid body exporter for the bin format (own file format).
 */
class RigidBodyExporter_BIN : public ExporterBase {
protected:
    bool m_isFirstFrame;
    std::string m_exportPath;

    void writeRigidBodies(const std::string& fileName);

public:
    explicit RigidBodyExporter_BIN(SimulatorBase* base);
    RigidBodyExporter_BIN(const RigidBodyExporter_BIN&) = delete;
    RigidBodyExporter_BIN& operator=(const RigidBodyExporter_BIN&) = delete;
    ~RigidBodyExporter_BIN() override;

    void init(const std::string& outputPath) override;
    void step(unsigned int frame) override;
    void reset() override;
    void setActive(const bool active) override;
};
}  // namespace vox