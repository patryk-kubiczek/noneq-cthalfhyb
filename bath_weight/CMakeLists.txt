cmake_minimum_required(VERSION 3.5)
project(bath_weight)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES BathWeight.cpp BathWeight.h main.cpp
                    ../model/cthyb/Model.cpp ../model/cthyb/Model.h
                    ../configuration/Configuration.cpp
                    ../random_number_generator/RandomNumberGenerator.cpp ../local_weight/LocalWeightWrapper.h ../model/ContourDeltaWrapper.h)

add_definitions(-DARMA_NO_DEBUG -DCTHYB_QMC)

add_executable(bath_weight ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(bath_weight ${TEST_LIBS})
