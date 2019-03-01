cmake_minimum_required(VERSION 3.5)
project(measurements)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES Measurements.h Measurements.cpp main.cpp
        measurements/DoubleOccupancy.h measurements/Occupancies.h measurements/InitialOccupancies.h measurements/MatsubaraGreenFunctions.h
        measurements/MeanComplexSign.h measurements/MeanOrder.h measurements/Occupancies.h measurements/OrderHistogram.h
        measurements/TimesStatistics.h measurements/EquilibriumGreenFunctions.h measurements/GreenFunctions.h
        measurements/SpinUp.h
        ../local_weight/cthalfhyb/LocalWeight.cpp ../model/cthalfhyb/Model.cpp
        ../discrete_bath/DiscreteBathWeight.cpp ../discrete_bath/DiscreteBath.cpp
        ../moves/Moves.cpp
        ../random_number_generator/RandomNumberGenerator.cpp ../configuration/Configuration.cpp ../bath_weight/BathWeight.cpp )

add_definitions(-DCTHALFHYB_QMC)
add_executable(measurements ${SOURCE_FILES})

set(TEST_LIBS ${LIBS} -lboost_iostreams -lboost_system -lboost_filesystem)

target_link_libraries(measurements ${TEST_LIBS})
