cmake_minimum_required(VERSION 3.15)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

project(reed-solomon LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories("include")

set(SOURCES "src")
set(LIBRS_SOURCES "src/librs")

add_library(rs STATIC
    "${LIBRS_SOURCES}/cyclotomic_coset.c"
    "${LIBRS_SOURCES}/fft.c"
    "${LIBRS_SOURCES}/gf65536.c"
    "${LIBRS_SOURCES}/reed_solomon.c"
    "${LIBRS_SOURCES}/symbol.c")

# ===== TESTING =====

# if(BUILD_TESTING)
#     enable_testing()
#     add_subdirectory(test)
# endif()