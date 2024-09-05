#!/bin/bash

# Number of CPUs
nnodes=6

## ----------------------------------------------
## Prepare time stamp
tstamp=$(date +'%y%m%d_%H%M%S')
echo TimeStamp: $tstamp

## ----------------------------------------------
## Copy input files to the working directory
cp ./input/input.dat ./inv/input.dat
cp ./input/drive.dat ./inv/drive.dat
cp ./input/station_info.dat ./inv/station_info.dat

## ----------------------------------------------
## Run PSI
if [ $nnodes = 1 ]; then
  # Run on single CPU
  nohup ./psi_sp > ${tstamp}.log 2>&1 &
else
  # Run on multiple CPU with MPI
  nohup mpirun -np $nnodes ./psi_pp > ${tstamp}.log 2>&1 &
fi

## ----------------------------------------------
## PID to file
echo PID: $!
echo $! > ${tstamp}.pid