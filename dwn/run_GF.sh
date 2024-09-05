#!/bin/bash

processors=8   #Number of parallel calculations (generally number of cores)

cat sources.dat | xargs -P$processors -n7 ./gr_nez
cat sources.dat | xargs -P$processors -n7 ./cnv_nez

./resort

echo "GF done"
