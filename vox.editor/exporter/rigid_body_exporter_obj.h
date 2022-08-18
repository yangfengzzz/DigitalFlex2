//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.editor/exporter/exporter_base.h"

namespace vox {
/** \brief Rigid body exporter for the OBJ format.
 */
class RigidBodyExporter_OBJ : public ExporterBase {
protected:
    bool m_isFirstFrame;
    std::string m_exportPath;

    void writeRigidBodies(unsigned int frame);

public:
    explicit RigidBodyExporter_OBJ(SimulatorBase* base);
    RigidBodyExporter_OBJ(const RigidBodyExporter_OBJ&) = delete;
    RigidBodyExporter_OBJ& operator=(const RigidBodyExporter_OBJ&) = delete;
    ~RigidBodyExporter_OBJ() override;

    void init(const std::string& outputPath) override;
    void step(unsigned int frame) override;
    void reset() override;
    void setActive(const bool active) override;
};
}  // namespace vox