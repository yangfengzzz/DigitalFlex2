cmake_minimum_required(VERSION 3.12)

project(apps.pbd LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.c
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/elastic_rods/position_based_elastic_rods_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/elastic_rods/position_based_elastic_rods_constraints.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/elastic_rods/position_based_elastic_rods_model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/elastic_rods/position_based_elastic_rods_tsc.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/fluid/fluid_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluid/fluid_model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluid/fluid_model_ts.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/generic_constraints/generic_constraints.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/generic_constraints/generic_constraints_model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/generic_constraints/generic_particle_constraints_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/generic_constraints/generic_rigid_body_constraints_app.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/stiff_rods/direct_position_based_solver_for_stiff_rods_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/stiff_rods/stretch_bending_twisting_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/stiff_rods/stiff_rods_scene_loader.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/bar_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/chain_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cloth_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cloth_collision_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cosserat_rods_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deformable_collision_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/joint_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_cloth_coupling_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_collision_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/scene_loader_app.cpp
        )


add_executable(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/AntTweakBar/include
        ../third_party/eigen
        ../third_party/glfw/include
        ../third_party/nlohmann/include
        ../third_party/spdlog/include
        ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(/System/Library/Frameworks)
find_library(APPKIT_LIBRARY AppKit)
find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} ${APPKIT_LIBRARY} ${OPENGL_LIBRARIES} glfw vox.base vox.pbd AntTweakBar OpenMP::OpenMP_CXX)

# Link third party libraries
#target_link_libraries(${PROJECT_NAME} PUBLIC
#        volk
#        ktx
#        stb
#        astc
#        imgui
#        glslang
#        SPIRV
#        vma
#        hwcpipe
#        spirv-cross-glsl
#        glslang-default-resource-limits
#        spdlog
#        ctpl
#        CLI11::CLI11)