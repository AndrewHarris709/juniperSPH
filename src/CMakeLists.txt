cmake_minimum_required(VERSION 3.16.3)
project(juniperSRC)

set(CMAKE_CXX_STANDARD 20)

add_library(juniperSRC STATIC kernel.cpp
        Simulation.cpp
        Simulation.h
        ParticleSet.cpp
        ParticleSet.h
        TreeNode.cpp
        TreeNode.h
        SimData.cpp
        SimData.h
        Integrator.h
        LeapfrogIntegrator.cpp
        LeapfrogIntegrator.h
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(juniperSRC PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(juniperSRC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
