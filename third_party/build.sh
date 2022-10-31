#!/bin/bash
set -x

# Jolt  ##################################
cd jolt/Build
mkdir build_release
cd build_release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j5
cd ..

mkdir build_debug
cd build_debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j5
cd ../../..

# GoogleTest ###############################
cd googletest
mkdir build
cd build
cmake ..
make -j5
cd ../..

set +x
