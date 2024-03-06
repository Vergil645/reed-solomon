#!/bin/bash

set -e

rm bin/*
rm lib/*
rm test/bin/*
rm test/lib/*

# Build
rm -rf build/debug/*
cmake -S. -Bbuild/debug -DBUILD_TESTING=true -DCMAKE_BUILD_TYPE=Debug
make -C build/debug