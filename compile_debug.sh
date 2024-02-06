#!/bin/bash

set -e

# Build
rm -rf build/*
cmake -S. -Bbuild -DBUILD_TESTING=true
make -C build