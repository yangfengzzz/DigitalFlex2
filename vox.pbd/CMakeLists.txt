cmake_minimum_required(VERSION 3.12)

project(vox.pbd LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${PROJECT_NAME} STATIC ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ../third_party/eigen
        ../third_party/nlohmann/include
        ../third_party/spdlog/include
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