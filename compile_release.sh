#!/bin/bash

rm bin/*
rm lib/*
rm -rf build/release/*

set -e

# Build
cmake -S. -Bbuild/release -DBUILD_TESTING=false -DCMAKE_BUILD_TYPE=Release -DADDITIONAL_C_FLAGS_RELEASE="$@"
make -C build/release