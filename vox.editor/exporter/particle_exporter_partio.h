//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <Partio.h>

#include <future>

#include "vox.editor/exporter/exporter_base.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Particle exporter for the partio format.
 */
class ParticleExporter_Partio : public ExporterBase {
protected:
    std::string m_exportPath;
    std::string m_particleFile;
    Partio::ParticlesDataMutable* m_particleData{};
    std::future<void> m_handle;

    void writeParticlesPartio(const std::string& fileName, FluidModel* model, unsigned int objId = 0xffffffff);

public:
    explicit ParticleExporter_Partio(SimulatorBase* base);
    ParticleExporter_Partio(const ParticleExporter_Partio&) = delete;
    ParticleExporter_Partio& operator=(const ParticleExporter_Partio&) = delete;
    ~ParticleExporter_Partio() override;

    void init(const std::string& outputPath) override;
    void step(unsigned int frame) override;
    void reset() override;
    void setActive(bool active) override;
};
}  // namespace vox