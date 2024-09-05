#!/bin/bash

## ----------------------------------------------
## Prepare time stamp
tstamp=$(date +'%y%m%d_%H%M%S')
echo TimeStamp: $tstamp

## ----------------------------------------------
## Copy input files to the working directory
cp ./input/afterprocess.dat ./inv/afterprocess.dat
FILE=var_Unezstd.tmp
if [ -f "$FILE" ]; then
    cp ${FILE} ./inv/${FILE}
fi

## ----------------------------------------------
## Run the afterprocessing of Trans-dimensional Inversion
nohup ./psi_ap > ${tstamp}.log 2>&1 &

## ----------------------------------------------
## PID to file
echo PID: $!
echo $! > ${tstamp}.pid