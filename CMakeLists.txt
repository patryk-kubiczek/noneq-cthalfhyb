cmake_minimum_required(VERSION 3.5)
project(noneq_physics)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost)
find_package(Armadillo)
add_definitions(-DARMA_DONT_USE_WRAPPER)
if (NOT ARMADILLO_FOUND)
    set(ARMADILLO_INCLUDE_DIRS $ENV{ARMADILLO_ROOT}/include)
endif (NOT ARMADILLO_FOUND)
find_package(BLAS)
find_package(LAPACK)
if (NOT LAPACK_FOUND OR NOT BLAS_FOUND)
    set(LAPACK_LIBRARIES "")
    set(BLAS_LIBRARIES "")
endif ()


message(STATUS "Compiler flags: " ${CMAKE_CXX_FLAGS})

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_subdirectory(equations_of_motion EXCLUDE_FROM_ALL)
add_subdirectory(input_generator EXCLUDE_FROM_ALL)

add_subdirectory(interpolation EXCLUDE_FROM_ALL)


add_subdirectory(model/cthalfhyb EXCLUDE_FROM_ALL)
add_subdirectory(model/cthyb EXCLUDE_FROM_ALL)


add_subdirectory(configuration EXCLUDE_FROM_ALL)
add_subdirectory(bath_weight EXCLUDE_FROM_ALL)
add_subdirectory(local_weight/cthyb EXCLUDE_FROM_ALL)
add_subdirectory(local_weight/cthalfhyb EXCLUDE_FROM_ALL)
add_subdirectory(moves/cthyb EXCLUDE_FROM_ALL)
add_subdirectory(moves/cthalfhyb EXCLUDE_FROM_ALL)
add_subdirectory(measurements EXCLUDE_FROM_ALL)
add_subdirectory(cthalfhyb_qmc EXCLUDE_FROM_ALL)
add_subdirectory(cthyb_qmc EXCLUDE_FROM_ALL)





add_subdirectory(dyson_equation EXCLUDE_FROM_ALL)
add_subdirectory(discrete_bath EXCLUDE_FROM_ALL)






