#!/bin/bash

rm bin/*
rm lib/*
rm test/bin/*
rm test/lib/*
rm -rf build/debug/*

set -e

# Build
echo
echo "========================"
echo "======= BUILDING ======="
echo "========================"
echo
cmake -S. -Bbuild/debug -DBUILD_TESTING=true -DCMAKE_BUILD_TYPE=Debug
make -C build/debug

# CppCheck
echo
echo "============================"
echo "===== RUNNING CppCheck ====="
echo "============================"
echo
cppcheck \
    --enable=all \
    --std=c11 \
    --inline-suppr \
    --inconclusive \
    --force \
    --suppress=missingIncludeSystem \
    -Iinclude -Itest/include \
    -isrc/rlc \
    src test/src