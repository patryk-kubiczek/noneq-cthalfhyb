cmake_minimum_required(VERSION 3.5)
project(model_cthalfhyb)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES Model.cpp Model.h main.cpp
        ../ContourDelta.h ../ContourPropagator.h ../analytic_expressions.h
        ../params.h ../contour_time.h)

add_executable(model_cthalfhyb ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(model_cthalfhyb ${TEST_LIBS})
