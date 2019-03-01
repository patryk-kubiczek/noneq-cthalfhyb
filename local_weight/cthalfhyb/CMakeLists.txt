cmake_minimum_required(VERSION 3.5)
project(local_weight_cthalfhyb)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES LocalWeight.h LocalWeight.cpp main.cpp
        ../../random_number_generator/RandomNumberGenerator.cpp ../../configuration/Configuration.cpp ../../model/cthalfhyb/Model.cpp ../../model/ContourPropagatorWrapper.h
        ../../discrete_bath/DiscreteBath.cpp ../../discrete_bath/DiscreteBathWeight.cpp)

add_definitions(-DCTHALFHYB_QMC)

add_executable(local_weight_cthalfhyb ${SOURCE_FILES})


set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(local_weight_cthalfhyb ${TEST_LIBS})
