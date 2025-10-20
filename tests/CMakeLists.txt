cmake_minimum_required(VERSION 3.16.3)
project(sph_tests)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Copy files into test directory
file(COPY ./files/test.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/files)
file(COPY ./files/hydro32_00020.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/files)
file(COPY ./files/kd_test.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/files)
file(COPY ./files/hydro32.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/files)
file(COPY ./files/test_integrate.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/files)

# Create a test for the kernel.cpp module:
add_executable(KernelTest test_kernel.cpp) # First create an initial executable from the test case code.
target_link_libraries(KernelTest PRIVATE gtest_main) # Link the test module to the gtest framework
target_link_libraries(KernelTest PRIVATE juniperSRC) # Add the kernel code modules to a `kernel_src` library
add_test(NAME Kernel_Tests COMMAND KernelTest) # Add the test as a new run configuration `KernelTest`

add_executable(IOTest test_io.cpp)
target_link_libraries(IOTest PRIVATE gtest_main gmock_main)
target_link_libraries(IOTest PRIVATE juniperSRC)
add_test(NAME IO_Tests COMMAND IOTest)

add_executable(NeighbourTest test_neighbours.cpp)
target_link_libraries(NeighbourTest PRIVATE gtest_main)
target_link_libraries(NeighbourTest PRIVATE juniperSRC)
add_test(NAME Neighbour_Tests COMMAND NeighbourTest)

add_executable(DensityTest test_density.cpp)
target_link_libraries(DensityTest PRIVATE gtest_main)
target_link_libraries(DensityTest PRIVATE juniperSRC)
add_test(NAME Density_Tests COMMAND DensityTest)
