cmake_minimum_required(VERSION 3.5)
project(moves_cthyb)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES ../Moves.h ../Moves.cpp ../main.cpp
        ../moves/ToggleFullEmptyLine.h ../moves/InsertSegment.h ../moves/RemoveSegment.h ../moves/ShiftSegmentEnd.h
        ../../bath_weight/BathWeight.cpp
        ../../local_weight/cthyb/LocalWeight.cpp ../../model/cthyb/Model.cpp
        ../../random_number_generator/RandomNumberGenerator.cpp ../../configuration/Configuration.cpp )

add_definitions(-DCTHYB_QMC)

add_executable(moves_cthyb ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(moves_cthyb ${TEST_LIBS})
