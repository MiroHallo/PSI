# --------------------------------------------------------------------------------
#
# Makefile for the Parametric Slip Inversion (PSI)
#
# Hallo, M., Gallovic, F. (2020): Bayesian self-adapting fault slip inversion 
#   with Green's functions uncertainty and application on the 2016 Mw7.1 Kumamoto
#   earthquake, Journal of Geophysical Research: Solid Earth, 125, e2019JB018703.
#
# --------------------------------------------------------------------------------
#
# Author: Miroslav Hallo
# Charles University, Faculty of Mathematics and Physics
# E-mail: hallo@karel.troja.mff.cuni.cz
# Revision 3/2019: The initial version
# 
# Copyright (C) 2019  Miroslav Hallo
# 
# This program is published under the GNU General Public License (GNU GPL).
# This program is free software: you can modify it and/or redistribute it
# or any derivative version under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
# This code is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
# and don't remove their names from the code.
# You should have received copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.
#
# --------------------------------------------------------------------------------
#
# Requirements: 
#    PT library, TIME_2D library, LAPACK or MKL, Fortran90 and MPI compilers
#
# Usage:
#    > make        # compile all binary files
#    > make clean  # remove all binary and module files
#
# --------------------------------------------------------------------------------

# SETTINGS (ETH)
# compiler without MPI
F90 = gfortran
# compiler with MPI
MF90 = mpif90
# flags
#FFLAGS = -O0 -cpp -Wall -pedantic -fbounds-check -fbacktrace -cpp
FFLAGS = -O3 -cpp
AFLAGS = -llapack -lblas
# path to the library folder
LP = ../lib/
# path to the parallel tempering library folder
PT = ../lib/pt/lib/

# SETTINGS (MFF)
#source load_impi
# compiler without MPI
#F90 = ifort
# compiler with MPI
#MF90 = mpiifort
# flags
#FFLAGS = -O3 -fpp
#AFLAGS = -mkl -DMKL
# path to the library folder
#LP = ../lib/
# path to the parallel tempering library folder
#PT = ../lib/pt/lib/

# SETTINGS (OVA)
#module load impi mkl
# compiler without MPI
#F90 = ifort
# compiler with MPI
#MF90 = mpiifort
# flags
#FFLAGS = -ipo -O3 -xCORE-AVX2 -fpp
#AFLAGS = -mkl -DMKL
# path to the library folder
#LP = ../lib/
# path to the parallel tempering library folder
#PT = ../lib/pt/lib/

# --------------------------------------------------------------------------------

all:: psi-stations psi-prepare psi-waves psi-sp psi-ap psi-pp
	@echo "Done."
	

# PSI STATIONS
psi-stations:
	$(F90) $(FFLAGS) -ostations stations.f90
	mv stations ../stations
	

# PSI PREPARE
psi-prepare:
	$(F90) $(FFLAGS) -oprepare prepare.f90
	mv prepare ../prepare
	

# PSI WAVES
psi-waves:
	$(F90) $(FFLAGS) -owaves $(LP)filters.for initInput.f90 waves.f90
	mv waves ../waves
	

# PSI SERIAL (without MPI)
psi-sp:
	$(F90) $(FFLAGS) -opsi_sp $(LP)filters.for $(LP)nr.for $(LP)time_2d.o $(PT)pt.serial.o invLib.f90 initInput.f90 initCore.f90 invExe.f90 initSynth.f90 psi.f90 $(AFLAGS)
	mv psi_sp ../psi_sp
	

# PSI RESULTS (without MPI)
psi-ap:
	$(F90) $(FFLAGS) -opsi_ap $(LP)filters.for $(LP)nr.for $(LP)time_2d.o $(PT)pt.serial.o invLib.f90 initInput.f90 initCore.f90 invExe.f90 initSynth.f90 finExe.f90 psiap.f90 $(AFLAGS)
	mv psi_ap ../psi_ap
	

# PSI PARALLEL (with MPI)
psi-pp:
	$(MF90) $(FFLAGS) -opsi_pp $(LP)filters.for $(LP)nr.for $(LP)time_2d.o $(PT)pt.mpi.o invLib.f90 initInput.f90 initCore.f90 invExe.f90 initSynth.f90 psi.f90 $(AFLAGS)
	mv psi_pp ../psi_pp
	

# Clean
clean: 
	-rm -f *.mod ../stations ../prepare ../waves ../psi_sp ../psi_ap ../psi_pp