cmake_minimum_required(VERSION 3.5)
project(local_weight_cthyb)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES LocalWeight.h LocalWeight.cpp main.cpp
        ../../random_number_generator/RandomNumberGenerator.cpp ../../configuration/Configuration.cpp ../../model/cthyb/Model.cpp ../../model/ContourPropagatorWrapper.h)

add_executable(local_weight_cthyb ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(local_weight_cthyb ${TEST_LIBS})
