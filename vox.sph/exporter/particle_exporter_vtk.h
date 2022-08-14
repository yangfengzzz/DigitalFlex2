//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <fstream>

#include "vox.sph/exporter/exporter_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Particle exporter for the VTK format.
 */
class ParticleExporter_VTK : public ExporterBase {
protected:
    std::string m_exportPath;
    std::ofstream* m_outfile;
    std::vector<std::string> m_attributes;

    void createParticleFile(const std::string& fileName, FluidModel* model);
    void writeParticles(const std::string& fileName, FluidModel* model, unsigned int objId = 0xffffffff);

    // VTK expects big endian
    template <typename T>
    inline void swapByteOrder(T* v) {
        constexpr size_t n = sizeof(T);
        auto* bytes = reinterpret_cast<uint8_t*>(v);
        for (unsigned int c = 0u; c < n / 2; c++) std::swap(bytes[c], bytes[n - c - 1]);
    }

public:
    explicit ParticleExporter_VTK(SimulatorBase* base);
    ParticleExporter_VTK(const ParticleExporter_VTK&) = delete;
    ParticleExporter_VTK& operator=(const ParticleExporter_VTK&) = delete;
    ~ParticleExporter_VTK() override;

    void init(const std::string& outputPath) override;
    void step(unsigned int frame) override;
    void reset() override;
    void setActive(bool active) override;
};
}  // namespace vox