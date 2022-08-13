cmake_minimum_required(VERSION 3.12)

project(vox.base LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utility/*.cpp)

add_library(${PROJECT_NAME} STATIC ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/eigen
        ../third_party/libmorton/include
        ../third_party/spdlog/include
        ${CMAKE_CURRENT_SOURCE_DIR})

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