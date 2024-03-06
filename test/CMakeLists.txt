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

set(RLC_TEST_SOURCES "src/rlc")
set(RS_TEST_SOURCES "src/rs")
set(LIBTESTUTIL_SOURCES "src/util")

add_library(testutil STATIC "${LIBTESTUTIL_SOURCES}/util.c")

# --- rs/gf65536

add_executable(test_rs_gf_mul_ee "${RS_TEST_SOURCES}/gf65536/test_gf_mul_ee.c")
target_link_libraries(test_rs_gf_mul_ee rs)

add_executable(test_rs_gf_div_ee "${RS_TEST_SOURCES}/gf65536/test_gf_div_ee.c")
target_link_libraries(test_rs_gf_div_ee rs)

# --- rs/cyclotomic_coset

add_executable(test_rs_cc_estimate_cosets_cnt "${RS_TEST_SOURCES}/cyclotomic_coset/test_cc_estimate_cosets_cnt.c")
target_link_libraries(test_rs_cc_estimate_cosets_cnt rs)

add_executable(test_rs_cc_select_cosets "${RS_TEST_SOURCES}/cyclotomic_coset/test_cc_select_cosets.c")
target_link_libraries(test_rs_cc_select_cosets rs)

add_executable(test_rs_cc_cosets_to_positions "${RS_TEST_SOURCES}/cyclotomic_coset/test_cc_cosets_to_positions.c")
target_link_libraries(test_rs_cc_cosets_to_positions rs testutil)

# --- rs/fft

# add_executable(test_rs_fft_transform "${RS_TEST_SOURCES}/fft/test_fft_transform.c")
# target_link_libraries(test_rs_fft_transform rs testutil)

# add_executable(test_rs_fft_partial_transform "${RS_TEST_SOURCES}/fft/test_fft_partial_transform.c")
# target_link_libraries(test_rs_fft_partial_transform rs testutil)

# --- rs/reed_solomon

# add_executable(test_rs_generate_repair_symbols "${RS_TEST_SOURCES}/reed_solomon/test_rs_generate_repair_symbols.c")
# target_link_libraries(test_rs_generate_repair_symbols rs testutil)

# add_executable(test_rs_restore_symbols "${RS_TEST_SOURCES}/reed_solomon/test_rs_restore_symbols.c")
# target_link_libraries(test_rs_restore_symbols rs testutil)

# --- rs

add_executable(test_rs_random_data "${RS_TEST_SOURCES}/test_random_data.c")
target_link_libraries(test_rs_random_data rs testutil)

# --- rlc

add_executable(test_rlc_random_data "${RLC_TEST_SOURCES}/test_random_data.c")
target_link_libraries(test_rlc_random_data rlc testutil)

# =========================
# ===== RUNNING TESTS =====
# =========================

# --- rs/gf65536

add_test(NAME test_rs_gf_mul_ee COMMAND test_rs_gf_mul_ee)
add_test(NAME test_rs_gf_div_ee COMMAND test_rs_gf_div_ee)

# --- rs/cyclotomic_coset

add_test(NAME test_rs_cc_estimate_cosets_cnt COMMAND test_rs_cc_estimate_cosets_cnt)
add_test(NAME test_rs_cc_select_cosets COMMAND test_rs_cc_select_cosets)
add_test(NAME test_rs_cc_cosets_to_positions COMMAND test_rs_cc_cosets_to_positions)

# --- rs/fft

# add_test(NAME test_rs_fft_transform COMMAND test_rs_fft_transform)
# add_test(NAME test_rs_fft_partial_transform COMMAND test_rs_fft_partial_transform)

# --- rs/reed_solomon

# add_test(NAME test_rs_generate_repair_symbols COMMAND test_rs_generate_repair_symbols)
# add_test(NAME test_rs_restore_symbols COMMAND test_rs_restore_symbols)

# --- rs

add_test(NAME test_rs_random_data COMMAND test_rs_random_data)

# --- rlc

add_test(NAME test_rlc_random_data COMMAND test_rlc_random_data)