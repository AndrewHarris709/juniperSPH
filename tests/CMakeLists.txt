cmake_minimum_required(VERSION 3.31)
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

# Create a test for the kernel.cpp module:
add_executable(KernelTest test_kernel.cpp
        test_neighbours.cpp) # First create an initial executable from the test case code.
target_link_libraries(KernelTest PRIVATE gtest_main) # Link the test module to the gtest framework
target_link_libraries(KernelTest PRIVATE juniperSRC) # Add the kernel code modules to a `kernel_src` library
add_test(NAME Kernel_Tests COMMAND KernelTest) # Add the test as a new run configuration `KernelTest`

# Create a test for the SimData.cpp module:
add_executable(ReaderTest test_reading.cpp)
target_link_libraries(ReaderTest PRIVATE gtest_main gmock_main)
target_link_libraries(ReaderTest PRIVATE juniperSRC)
add_test(NAME Reader_Tests COMMAND ReaderTest)

add_executable(NeighbourTest test_neighbours.cpp)
target_link_libraries(NeighbourTest PRIVATE gtest_main)
target_link_libraries(NeighbourTest PRIVATE juniperSRC)
add_test(NAME Neighbour_Tests COMMAND NeighbourTest)

add_executable(DensityTest test_density.cpp)
target_link_libraries(DensityTest PRIVATE gtest_main)
target_link_libraries(DensityTest PRIVATE juniperSRC)
add_test(NAME Density_Tests COMMAND DensityTest)
