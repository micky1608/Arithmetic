cmake_minimum_required(VERSION 3.12)
project(Arithmetic C)

set(CMAKE_C_STANDARD 11)

add_executable(Arithmetic main.c arithmetic.c pol_bigint_arithmetic.c pol_bigint_arithmetic.h Util.c Util.h pol_bigfloat_arithmetic.c pol_bigfloat_arithmetic.h Euclidean.c Euclidean.h pol_bignumber.c pol_bignumber.h)
target_link_libraries(Arithmetic gmp m)