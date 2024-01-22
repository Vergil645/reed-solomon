set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

project(reed-solomon_test LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCES "src")
set(LIBRS_SOURCES "src/librs")

add_executable(test_gf_mul_ee "${LIBRS_SOURCES}/gf65536/test_gf_mul_ee.c")
target_link_libraries(test_gf_mul_ee rs)

add_executable(test_gf_div_ee "${LIBRS_SOURCES}/gf65536/test_gf_div_ee.c")
target_link_libraries(test_gf_div_ee rs)

# ===== TESTS =====

add_test(NAME test_gf_mul_ee COMMAND test_gf_mul_ee)
add_test(NAME test_gf_div_ee COMMAND test_gf_div_ee)