#!/usr/bin/env bash

export BUILD_TYPE=Release
export CXXFLAGS="-DNDEBUG -O3 -std=gnu++14 -fext-numeric-literals"
export LINK_LIBRARIES=ON
export CC=mpicc
export CXX=mpicxx

cmake -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_CXX_FLAGS_RELEASE=$CXXFLAGS \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DBUILD_SHARED_LIBS=$LINK_LIBRARIES .
make eom -j4
make input_generator -j4

make cthyb -j4
make cthalfhyb -j4

make cthyb_qmc -j4
make cthalfhyb_qmc -j4

pushd ./python
python3 setup_input.py build_ext --inplace

QMC=CTHALFHYB python3 setup_qmc.py build_ext --inplace
QMC=CTHYB python3 setup_qmc.py build_ext --inplace

popd
