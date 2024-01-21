#!/bin/bash

rm -rf build && \
    cmake -S. -Bbuild && \
    make -C build