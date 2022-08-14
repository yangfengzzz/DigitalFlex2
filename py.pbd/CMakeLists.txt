cmake_minimum_required(VERSION 3.12)

project(py.pbd LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

pybind11_add_module(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/eigen
        ../third_party/nlohmann/include
        ../third_party/spdlog/include
        ../third_party/pybind11/include
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