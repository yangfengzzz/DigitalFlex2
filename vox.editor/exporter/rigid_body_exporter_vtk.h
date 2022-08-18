//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.editor/exporter/exporter_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Rigid body exporter for the VTK format.
 */
class RigidBodyExporter_VTK : public ExporterBase {
protected:
    bool m_isFirstFrame;
    std::string m_exportPath;

    void writeRigidBodies(unsigned int frame);

    // VTK expects big endian
    template <typename T>
    inline void swapByteOrder(T* v) {
        constexpr size_t n = sizeof(T);
        auto* bytes = reinterpret_cast<uint8_t*>(v);
        for (unsigned int c = 0u; c < n / 2; c++) std::swap(bytes[c], bytes[n - c - 1]);
    }

public:
    explicit RigidBodyExporter_VTK(SimulatorBase* base);
    RigidBodyExporter_VTK(const RigidBodyExporter_VTK&) = delete;
    RigidBodyExporter_VTK& operator=(const RigidBodyExporter_VTK&) = delete;
    ~RigidBodyExporter_VTK() override;

    void init(const std::string& outputPath) override;
    void step(unsigned int frame) override;
    void reset() override;
    void setActive(const bool active) override;
};
}  // namespace vox