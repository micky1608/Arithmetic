cmake_minimum_required(VERSION 3.12)
project(Arithmetic C)

set(CMAKE_C_STANDARD 11)

add_executable(Arithmetic main.c arithmetic.c polynomial_arithmetic.c polynomial_arithmetic.h Util.c Util.h)
target_link_libraries(Arithmetic gmp)