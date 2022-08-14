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

#include "vox.base/common.h"
#include "vox.pbd/simulation_model.h"

// PYBIND11_MAKE_OPAQUE(std::vector<Vector3r>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::TriangleModel*>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::TetModel*>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::Constraint*>)
PYBIND11_MAKE_OPAQUE(std::vector<vox::RigidBody*>)