cmake_minimum_required(VERSION 3.31)
project(juniperEXE)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(juniperEXE main.cpp)
target_link_libraries(juniperEXE PUBLIC juniperSRC)
