cmake_minimum_required(VERSION 3.12)
project(Arithmetic C)

set(CMAKE_C_STANDARD 11)

add_executable(Arithmetic main.c arithmetic.c
        pol_bigint_arithmetic.c pol_bigint_arithmetic.h
        Util.c Util.h pol.c pol.h
        pol_bigfloat_arithmetic.c pol_bigfloat_arithmetic.h
        Euclidean.c Euclidean.h matrix.h matrix.c
        pol_bignumber.c pol_bignumber.h
        matrix_bigint_arithmetic.c matrix_bigint_arithmetic.h
        matrix_bigQ_arithmetic.c matrix_bigQ_arithmetic.h CRT_bigint.c CRT_bigint.h Lagrange.c Lagrange.h
        pol_bigQ_arithmetic.c pol_bigQ_arithmetic.h matrix_double.c matrix_double.h
        Complex.c Complex.h FFT.c FFT.h test.c test.h )

target_link_libraries(Arithmetic gmp mpfr m mpc)