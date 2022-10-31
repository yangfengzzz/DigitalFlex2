#!/bin/bash
set -x

# Jolt  #####################################
cd jolt/Build
rm -rf build_release
rm -rf build_debug
cd ../..

# GoogleTest ##############################
cd googletest
rm -rf build
cd ..


set +x
