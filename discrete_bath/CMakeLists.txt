cmake_minimum_required(VERSION 3.5)
project(discrete_bath)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES DiscreteBath.cpp DiscreteBath.h main.cpp
        DiscreteBathWeight.h DiscreteBathWeight.cpp
        ../configuration/Configuration.cpp ../random_number_generator/RandomNumberGenerator.cpp)

add_definitions(-DARMA_NO_DEBUG -DCTHALFHYB_QMC)

add_executable(discrete_bath ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)


target_link_libraries(discrete_bath ${TEST_LIBS})
