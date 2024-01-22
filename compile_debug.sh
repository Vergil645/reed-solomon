#!/bin/bash

rm -rf build && \
    cmake -S. -Bbuild -DBUILD_TESTING=true && \
    make -C build && \
    ctest --test-dir build