cmake_minimum_required(VERSION 3.5)
project(equations_of_motion)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES
        main.cpp
        EquationsOfMotion.h EquationsOfMotion.cpp
        Solver.h
        time_dependent_hamiltonian.cpp time_dependent_hamiltonian.h)

add_library(eom ${SOURCE_FILES})
add_executable(equations_of_motion ${SOURCE_FILES})

target_link_libraries(eom ${LIBS})
target_link_libraries(equations_of_motion ${LIBS})