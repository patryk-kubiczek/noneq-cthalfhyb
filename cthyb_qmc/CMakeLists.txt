cmake_minimum_required(VERSION 3.5)
project(cthyb_qmc)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES
		
        ../model/cthyb/Model.h ../model/cthyb/Model.cpp
		../model/params.h 
        ../model/contour_time.h 
		../model/ContourDelta.h 
        ../model/ContourPropagator.h
		../model/analytic_expressions.h
        ../model/ContourDeltaWrapper.h
        ../model/ContourPropagatorWrapper.h
        
        ../configuration/Configuration.h ../configuration/Configuration.cpp

        ../bath_weight/BathWeight.h ../bath_weight/BathWeight.cpp

        ../local_weight/LocalWeightBase.h
        ../local_weight/LocalWeightWrapper.h
        
        ../local_weight/cthyb/LocalWeight.h ../local_weight/cthyb/LocalWeight.cpp

        ../moves/Moves.h ../moves/Moves.cpp
        ../moves/moves/ToggleFullEmptyLine.h
        ../moves/moves/InsertSegment.h
        ../moves/moves/RemoveSegment.h
        ../moves/moves/ShiftSegmentEnd.h

        ../measurements/Measurements.h ../measurements/Measurements.cpp
        ../measurements/measurements/MeanComplexSign.h
        ../measurements/measurements/MeanOrder.h
        ../measurements/measurements/OrderHistogram.h
        ../measurements/measurements/TimesStatistics.h
        ../measurements/measurements/InitialOccupancies.h
        ../measurements/measurements/Occupancies.h
        ../measurements/measurements/GreenFunctions.h
        ../measurements/measurements/EquilibriumGreenFunctions.h
        ../measurements/measurements/MatsubaraGreenFunctions.h
        ../measurements/measurements/SpinUp.h

        ../random_number_generator/RandomNumberGenerator.h ../random_number_generator/RandomNumberGenerator.cpp

        ../interpolation/Interpolation.h
		../data_types/Tetracube.h

        ../auxiliary_functions/get_permutation_sign.h
        ../auxiliary_functions/container_functions.h
        ../auxiliary_functions/matrix_functions.h
        ../auxiliary_functions/constants.h
        ../auxiliary_functions/complex_numbers.h
        ../auxiliary_functions/element_wise_functions.h
        ../auxiliary_functions/mpi_wrappers.h
        ../auxiliary_functions/double_to_int.h
		../auxiliary_functions/save_functions.h )

#add_definitions(-DARMA_DONT_USE_WRAPPER)
#set(ENV{NVBLAS_CONFIG_FILE} "${dir}/bin/nvblas.conf")
#set(LIBS ${LIBS} -L/usr/local/cuda-8.0/lib64/libnvblas.so)

add_definitions(-DCTHYB_QMC)

add_executable(cthyb_qmc QMC.h QMC.cpp main.cpp ${SOURCE_FILES})
target_link_libraries(cthyb_qmc ${LIBS})
target_compile_definitions(cthyb_qmc PUBLIC)

add_library(cthyb ${SOURCE_FILES})
target_link_libraries(cthyb ${LIBS})
target_compile_definitions(cthyb PUBLIC -DUSE_MPI)

# Run the executable with LD_PRELOAD:
# LD_PRELOAD=/usr/local/cuda-8.0/lib64/libnvblas.so ./noneq_qmc
