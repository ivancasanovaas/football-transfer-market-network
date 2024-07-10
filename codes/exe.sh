#!/bin/bash

echo
gfortran -c structural_properties.f95
echo
gfortran -c equilibrium_models.f95
echo
gfortran -c SIS_dynamics.f95
echo
gfortran -c main.f95
echo
gfortran -o main main.o structural_properties.o equilibrium_models.o SIS_dynamics.o
echo
./main ../network/edges.txt ../results/dat
