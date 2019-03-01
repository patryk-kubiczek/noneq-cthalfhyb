cmake_minimum_required(VERSION 3.5)
project(configuration)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES Configuration.cpp Configuration.h main.cpp
                    ../random_number_generator/RandomNumberGenerator.cpp)

add_definitions(-DARMA_NO_DEBUG)

add_executable(configuration ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)


target_link_libraries(configuration ${TEST_LIBS})
