
cmake_minimum_required(VERSION 2.6)

project(s2mini)

set(s2mini_version_major 2)
set(s2mini_version_minor 0)

FIND_PACKAGE(MPI REQUIRED)
INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER(${MPI_COMPILER} "MPI C++ Compiler")

# OpenMP
FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)

    # hard-wired -mp for PGI compiler
    set (OpenMP_C_FLAGS "-mp")
    set (OpenMP_CXX_FLAGS "-mp")

    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()


set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

#OpenSSL Libraries
#find_package(OpenSSL REQUIRED)

# for vim ycm plugin
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

include_directories(
    ${PROJECT_SOURCE_DIR}/install/include
    ${PROJECT_SOURCE_DIR}/src
)

link_directories(
    ${PROJECT_SOURCE_DIR}/install/lib
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)
add_subdirectory(examples)