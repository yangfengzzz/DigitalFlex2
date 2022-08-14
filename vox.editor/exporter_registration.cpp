//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.editor/exporter/exporter_base.h"
#include "vox.editor/exporter/particle_exporter_partio.h"
#include "vox.editor/exporter/particle_exporter_vtk.h"
#include "vox.editor/exporter/rigid_body_exporter_bin.h"
#include "vox.editor/exporter/rigid_body_exporter_obj.h"
#include "vox.editor/exporter/rigid_body_exporter_vtk.h"

using namespace vox;

void SimulatorBase::createExporters() {
    addParticleExporter("enablePartioExport", "Partio Exporter", "Enable/disable partio export.",
                        new ParticleExporter_Partio(this));
    addParticleExporter("enableVTKExport", "VTK Exporter", "Enable/disable VTK export.",
                        new ParticleExporter_VTK(this));

    addRigidBodyExporter("enableRigidBodyExport", "Rigid Body Exporter", "Enable/disable rigid body BIN export.",
                         new RigidBodyExporter_BIN(this));
    addRigidBodyExporter("enableRigidBodyOBJExport", "Rigid Body OBJ Exporter", "Enable/disable rigid body OBJ export.",
                         new RigidBodyExporter_OBJ(this));
    addRigidBodyExporter("enableRigidBodyVTKExport", "Rigid Body VTK Exporter", "Enable/disable rigid body VTK export.",
                         new RigidBodyExporter_VTK(this));
}