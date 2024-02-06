set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

project(reed-solomon_test LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ==========================
# ===== BUILDING TESTS =====
# ==========================

include_directories("include")

set(SOURCES "src")
set(LIBRS_TEST_SOURCES "src/librs")

add_library(util STATIC "${SOURCES}/util.c")

# --- librs/gf65536

add_executable(test_gf_mul_ee "${LIBRS_TEST_SOURCES}/gf65536/test_gf_mul_ee.c")
target_link_libraries(test_gf_mul_ee rs)

add_executable(test_gf_div_ee "${LIBRS_TEST_SOURCES}/gf65536/test_gf_div_ee.c")
target_link_libraries(test_gf_div_ee rs)

# --- librs/cyclotomic_coset

add_executable(test__cc_get_cosets_cnt "${LIBRS_TEST_SOURCES}/cyclotomic_coset/test__cc_get_cosets_cnt.c")
target_link_libraries(test__cc_get_cosets_cnt rs)

add_executable(test_cc_estimate_cosets_cnt "${LIBRS_TEST_SOURCES}/cyclotomic_coset/test_cc_estimate_cosets_cnt.c")
target_link_libraries(test_cc_estimate_cosets_cnt rs)

add_executable(test_cc_select_cosets "${LIBRS_TEST_SOURCES}/cyclotomic_coset/test_cc_select_cosets.c")
target_link_libraries(test_cc_select_cosets rs)

add_executable(test_cc_cosets_to_positions "${LIBRS_TEST_SOURCES}/cyclotomic_coset/test_cc_cosets_to_positions.c")
target_link_libraries(test_cc_cosets_to_positions rs util)

# --- librs/fft

add_executable(test_fft_transform "${LIBRS_TEST_SOURCES}/fft/test_fft_transform.c")
target_link_libraries(test_fft_transform rs util)

add_executable(test_fft_partial_transform "${LIBRS_TEST_SOURCES}/fft/test_fft_partial_transform.c")
target_link_libraries(test_fft_partial_transform rs util)

# --- librs/reed_solomon

add_executable(test_rs_generate_repair_symbols "${LIBRS_TEST_SOURCES}/reed_solomon/test_rs_generate_repair_symbols.c")
target_link_libraries(test_rs_generate_repair_symbols rs util)

# =========================
# ===== RUNNING TESTS =====
# =========================

# --- librs/gf65536

add_test(NAME test_gf_mul_ee COMMAND test_gf_mul_ee)
add_test(NAME test_gf_div_ee COMMAND test_gf_div_ee)

# --- librs/cyclotomic_coset

add_test(NAME test__cc_get_cosets_cnt COMMAND test__cc_get_cosets_cnt)
add_test(NAME test_cc_estimate_cosets_cnt COMMAND test_cc_estimate_cosets_cnt)
add_test(NAME test_cc_select_cosets COMMAND test_cc_select_cosets)
add_test(NAME test_cc_cosets_to_positions COMMAND test_cc_cosets_to_positions)

# --- librs/fft

add_test(NAME test_fft_transform COMMAND test_fft_transform)
add_test(NAME test_fft_partial_transform COMMAND test_fft_partial_transform)

# --- librs/reed_solomon

add_test(NAME test_rs_generate_repair_symbols COMMAND test_rs_generate_repair_symbols)