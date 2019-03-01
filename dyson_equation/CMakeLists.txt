cmake_minimum_required(VERSION 3.5)
project(dyson_equation)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES DysonEquation.cpp DysonEquation.h main.cpp)


add_executable(dyson_eq ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(dyson_eq ${TEST_LIBS})
