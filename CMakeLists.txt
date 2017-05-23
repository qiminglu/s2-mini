
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
    set (OpenMP_PGI_C_FLAGS "-mp")
    set (OpenMP_PGI_CXX_FLAGS "-mp")

    set (OpenMP_GCC_C_FLAGS "-fopenmp")
    set (OpenMP_GCC_CXX_FLAGS "-fopenmp")

endif ()

set(OpenACC_PGI_C_FLAGS   "-acc -ta=tesla:cc60,time -Minfo=accel")
set(OpenACC_PGI_CXX_FLAGS "-acc -ta=tesla:cc60,time -Minfo=accel")

set(OpenACC_GCC_GXX_FLAGS "-fopenacc -foffload=nvptx-none ")


#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS "-std=c++11 -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${OpenACC_PGI_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS "-std=c++11 -fopenacc -foffload=nvptx-none -foffload=\"-O3\" -O3")
#set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -O3")

#set(CMAKE_CXX_FLAGS "-std=c++11 -acc -ta=tesla:cc60,time -Minfo=accel -O3")
set(CMAKE_CXX_FLAGS "-std=c++11 -mp -O3")

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
