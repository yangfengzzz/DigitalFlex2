cmake_minimum_required(VERSION 3.12)

project(vox.sph LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/drag/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/elasticity/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/surface_tension/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/viscosity/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/vorticity/*.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/dfsph/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/icsph/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/iisph/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pcisph/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/wcsph/*.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/pbf/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/projective_fluids/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sampling/*.cpp)

add_library(${PROJECT_NAME} STATIC ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/eigen
        ../third_party/nlohmann/include
        ../third_party/spdlog/include
        ../third_party/tinyexpr
        ../third_party/partio/src/lib
        ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

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