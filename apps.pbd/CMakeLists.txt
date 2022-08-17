cmake_minimum_required(VERSION 3.12)

project(apps.pbd LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.c
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/AntTweakBar/include
        ../third_party/eigen
        ../third_party/glfw/include
        ../third_party/nlohmann/include
        ../third_party/spdlog/include
        ${CMAKE_CURRENT_SOURCE_DIR})

include_directories ( /System/Library/Frameworks )
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