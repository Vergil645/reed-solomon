#!/bin/bash

set -e


perf record -a -g -- ./bin/compare_codes
perf script > out.perf

/home/matveyk/software/flamegraph/stackcollapse-perf.pl out.perf > out.folded
/home/matveyk/software/flamegraph/flamegraph.pl "$@" out.folded > out.svg