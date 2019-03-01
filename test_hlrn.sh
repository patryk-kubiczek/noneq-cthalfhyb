#!/bin/bash

#PBS -N noneq_qmc_test
#PBS -j oe
#PBS -l walltime=0:05:00
#PBS -l nodes=1:ppn=24
#PBS -l feature=mpp:test   		

mkdir -p input
mkdir -p output
mkdir -p output/CTHYB
mkdir -p output/CTHALFHYB
mkdir -p plots

QMC=CTHYB GENERATE_INPUT=1 aprun -b -n 24 python3 test.py > test_cthyb.out
QMC=CTHALFHYB GENERATE_INPUT=0 aprun -b -n 24 python3 test.py > test_cthalfhyb.out

python3 analyze_data.py

