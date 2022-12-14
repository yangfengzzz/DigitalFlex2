# Copyright (c) 2020-2021, Arm Limited and Contributors
#
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 the "License";
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

cmake_minimum_required(VERSION 3.12)

project(Flex-cpp)

# create output folder
file(MAKE_DIRECTORY output)

# Add path for local cmake scripts
list(APPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/bldsys
    ${CMAKE_CURRENT_SOURCE_DIR}/bldsys/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/bldsys/cmake/module)

# Precompiled header macro. Parameters are source file list and filename for pch cpp file.
macro(vulkan_samples_pch SRCS PCHCPP)
    if(MSVC AND CMAKE_GENERATOR MATCHES "^Visual Studio")
        set(PCH_NAME "$(IntDir)\\pch.pch")
        # make source files use/depend on PCH_NAME
        set_source_files_properties(${${SRCS}} PROPERTIES COMPILE_FLAGS "/Yupch.h /FIpch.h /Fp${PCH_NAME} /Zm300" OBJECT_DEPENDS "${PCH_NAME}")
        # make PCHCPP file compile and generate PCH_NAME
        set_source_files_properties(${PCHCPP} PROPERTIES COMPILE_FLAGS "/Ycpch.h /Fp${PCH_NAME} /Zm300" OBJECT_OUTPUTS "${PCH_NAME}")
        list(APPEND ${SRCS} "${PCHCPP}")
    endif()
endmacro(vulkan_samples_pch)

include(utils)
include(global_options)
include(sample_helper)
include(check_atomic)

find_package(pybind11)
include(pybind11Common)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif ()

####################################################
# Add third party libraries
add_subdirectory(third_party)
add_subdirectory(third_party/AntTweakBar)

# Add base library
add_subdirectory(vox.base)
add_subdirectory(test.base)

## Add pbd library
add_subdirectory(vox.pbd)
add_subdirectory(py.pbd)
add_subdirectory(apps.pbd)

# Add sph library
add_subdirectory(vox.sph)
add_subdirectory(py.sph)
add_subdirectory(test.sph)
#add_subdirectory(vox.editor)