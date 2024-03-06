#!/bin/bash

set -e

if [[ -d external/flec-moepgf ]]; then
    pushd external/flec-moepgf
    autoreconf -fi
    ./configure
    make
    cp .libs/libmoepgf.a ../
    popd
fi