cmake_minimum_required(VERSION 3.22)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

project(reed-solomon LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ================================
# ===== BUILDING THE PROJECT =====
# ================================

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fanalyzer")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${ADDITIONAL_C_FLAGS_RELEASE} -Wall -DNDEBUG")

include_directories("include")

set(EXTERNAL "external")
set(SOURCES "src")
set(LIBMEMORY_SOURCES "src/memory")
set(LIBRLC_SOURCES "src/rlc")
set(LIBRS_SOURCES "src/rs")

add_library(memory STATIC
    "${LIBMEMORY_SOURCES}/seq.c"
    "${LIBMEMORY_SOURCES}/symbol.c")

add_library(rlc STATIC
    "${LIBRLC_SOURCES}/equation.c"
    "${LIBRLC_SOURCES}/gf256.c"
    "${LIBRLC_SOURCES}/rlc.c"
    "${LIBRLC_SOURCES}/system.c")
target_link_libraries(rlc memory)

add_library(rs STATIC
    "${LIBRS_SOURCES}/cyclotomic_coset.c"
    "${LIBRS_SOURCES}/fft.c"
    "${LIBRS_SOURCES}/gf65536.c"
    "${LIBRS_SOURCES}/reed_solomon.c")
target_link_libraries(rs memory)

add_executable(compare_op_gf256 "${SOURCES}/compare_op_gf256.c")
target_link_libraries(compare_op_gf256 rlc)

add_executable(compare_codes "${SOURCES}/compare_codes.c")
target_link_libraries(compare_codes rs rlc)

add_executable(run_enc_dec "${SOURCES}/run_enc_dec.c")
target_link_libraries(run_enc_dec rs rlc)

# ===================
# ===== TESTING =====
# ===================

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()