cmake_minimum_required(VERSION 3.1)

project(nbody)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3 -Wall")

set(HEADERS lib/nBodySim.hpp
            lib/node.hpp
            lib/tree.hpp
            lib/utils.hpp
            )

set(SOURCES lib/nBodySim.cpp
            lib/node.cpp
            lib/tree.cpp
            lib/utils.cpp
            )

add_library(nbodyLib ${HEADERS} ${SOURCES})
set_target_properties(nbodyLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(nbodyLib PUBLIC lib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nbodyLib PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(nbody main.cpp)
target_link_libraries(nbody nbodyLib)