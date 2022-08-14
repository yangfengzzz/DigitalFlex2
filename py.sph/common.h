//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "vox.editor/position_based_dynamics_wrapper/pbd_wrapper.h"
#include "vox.sph/animation_field.h"
#include "vox.sph/emitter.h"
#include "vox.sph/obj_loader.h"
#include "vox.sph/utilities/scene_loader.h"

PYBIND11_MAKE_OPAQUE(std::vector<vox::Emitter *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::AnimationField *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::FieldDescription>)
PYBIND11_MAKE_OPAQUE(std::vector<std::array<float, 3>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::array<float, 2>>)
PYBIND11_MAKE_OPAQUE(std::vector<Vector3r>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::MeshFaceIndices>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::BoundaryData *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::FluidData *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::FluidBlock *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::EmitterData *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::AnimationFieldData *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::utility::SceneLoader::MaterialData *>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::PBDWrapper::RBData>)