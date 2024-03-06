#!/bin/bash

set -e

rm bin/*
rm lib/*

# Build
rm -rf build/release/*
cmake -S. -Bbuild/release -DBUILD_TESTING=false -DCMAKE_BUILD_TYPE=Release -DADDITIONAL_C_FLAGS_RELEASE="$@"
make -C build/release