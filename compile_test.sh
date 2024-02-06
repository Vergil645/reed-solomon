#!/bin/bash

set -e

# Build
./compile_debug.sh

# Run Ctest
echo
echo "========================="
echo "===== RUNNING CTEST ====="
echo "========================="
echo
ctest --test-dir build --stop-on-failure

# Run tests with Valgind

echo
echo "============================"
echo "===== RUNNING VALGRIND ====="
echo "============================"
echo
find test/bin -name "test_*" | xargs -t -l1 valgrind -s -q --leak-check=full