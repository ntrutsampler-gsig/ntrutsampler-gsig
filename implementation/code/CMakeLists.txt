cmake_minimum_required(VERSION 3.11)

project(sig)

# Use C11 and crash if that is not available/possible
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(test)
add_executable(bench)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")
set(BENCH_DIR "${CMAKE_SOURCE_DIR}/bench")
set(INCLUDE_DIR
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/include/arith")

set(SOURCES_SIG
  "${SRC_DIR}/arith/arith.c"
  "${SRC_DIR}/arith/poly_q.c"
  "${SRC_DIR}/arith/poly_q_sampling.c"
  "${SRC_DIR}/arith/poly_real.c"
  "${SRC_DIR}/precomputations.c"
  "${SRC_DIR}/fips202.c"
  "${SRC_DIR}/fpr.c"
  "${SRC_DIR}/fft.c"
  "${SRC_DIR}/random.c"
  "${SRC_DIR}/randombytes.c"
  "${SRC_DIR}/sampling.c"
  "${SRC_DIR}/sign.c")

target_sources(test PUBLIC
  "${TEST_DIR}/main.c"
  "${SOURCES_SIG}")

target_sources(bench PUBLIC
  "${BENCH_DIR}/main.c"
  "${BENCH_DIR}/benchmark.c"
  "${BENCH_DIR}/bench_sign.c"
  "${BENCH_DIR}/bench_sampler.c"
  "${SOURCES_SIG}")

target_include_directories(test PUBLIC ${INCLUDE_DIR})
target_include_directories(bench PUBLIC ${INCLUDE_DIR} ${BENCH_DIR})

set(WARN_ERROR_FLAGS
  -Wall -Wextra -Wundef -Wshadow -Wpointer-arith
  -Wcast-align -Wstrict-prototypes -Wwrite-strings -Waggregate-return
  -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code)

set(SANITIZER_FLAGS
  -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer
  -fsanitize=undefined -fsanitize=float-divide-by-zero
  -fsanitize=float-cast-overflow)

set(TEST_OPTS
  -march=native
  -g
  -O3
  )

set(BENCH_OPTS
  -march=native
  -O3
  )

target_compile_options(test PUBLIC ${WARN_ERROR_FLAGS} ${TEST_OPTS} ${SANITIZER_FLAGS})
target_link_options(test PUBLIC ${TEST_OPTS} ${SANITIZER_FLAGS})

target_compile_definitions(bench PUBLIC REMOVE_ASSERTS=1)
target_compile_options(bench PUBLIC ${WARN_ERROR_FLAGS} ${BENCH_OPTS})
target_link_options(bench PUBLIC ${BENCH_OPTS})

target_link_libraries(test PUBLIC "${CMAKE_SOURCE_DIR}/libs/flint-3.5.0/libflint.so")
target_link_libraries(test PUBLIC m)

target_link_libraries(bench PUBLIC "${CMAKE_SOURCE_DIR}/libs/flint-3.5.0/libflint.so")
target_link_libraries(bench PUBLIC m)
