cmake_minimum_required(VERSION 3.5)
project(input_generator)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)

set(EOM_DIR ${PARENT_DIR}/equations_of_motion)
if (BUILD_SHARED_LIBS)
    set(EOM_LIB ${PARENT_DIR}/equations_of_motion/lib/libeom.so)
else()
    set(EOM_LIB ${PARENT_DIR}/equations_of_motion/lib/libeom.a)
endif ()

set(SOURCE_FILES InputGenerator.cpp InputGenerator.h)

add_executable(input_generator main.cpp ${SOURCE_FILES})
target_link_libraries(input_generator ${EOM_LIB} ${LIBS})

