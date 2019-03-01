#!/usr/bin/env bash

mkdir -p input
mkdir -p output
mkdir -p output/CTHYB
mkdir -p output/CTHALFHYB
mkdir -p plots

QMC=CTHYB GENERATE_INPUT=1 mpirun -np 4 python3 test.py
QMC=CTHALFHYB GENERATE_INPUT=0 mpirun -np 4 python3 test.py

python3 analyze_data.py
