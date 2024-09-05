# Parametric Slip Inversion (PSI)
***************************************

  Open-source functions for the Trans-dimensional Bayesian fault slip
  inversion folowing Hallo and Gallovic (2020).

1 METHODOLOGY
===================

  Hallo, M., Gallovic, F. (2020). Bayesian self-adapting fault slip
inversion with Green's functions uncertainty and application on the
2016 Mw7.1 Kumamoto earthquake, Journal of Geophysical Research: 
Solid Earth, 125, e2019JB018703.

2 COPYRIGHT
===================

Copyright (C) 2017-2019  Miroslav Hallo

This program is published under the GNU General Public License (GNU GPL).

This program is free software: you can modify it and/or redistribute it
or any derivative version under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
and don't remove their names from the code.

You should have received copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

3 REQUIREMENTS
===================

  a) Compiled "pt.f90" 
    - Parallel Tempering library (M.Sambridge)
    - http://www.iearth.org.au/codes/ParallelTempering/
  b) Compiled "time_2d.c" 
    - Finite-differences computation of 2D travel time library (P.Podvin)
    - https://github.com/fgallovic/RIKsrf/blob/master/src-RIKsrf/Time_2d.c
  c) Compiled "gr_nez.for" and "cnv_nez.for"
    - Discrete wavenumber method for Green's functions computation (M.Bouchon)
    - https://github.com/fgallovic/LinSlipInv/tree/master/src-dwn
  d) LINUX/UNIX machine with LAPACK or MKL
  e) Fortran90 and MPI compilers
  f) MATLAB for ploting results

4 PACKAGE STRUCTURE
===================

  <dwn>       - Directory containing Green's functions computation
  <examples>  - Directory containing examples of input files
  <input>     - Directory with input files
  <inv>       - Work directory for the ongoing inversion
  <lib>       - Directory with compiled libraries and additional functions
  <src>       - Directory with source codes of the PSI
  plot_prepare.m            - Plots prepared suf-faults for Green's functions
  plot_psi_fit.m            - Plots final data fit
  plot_psi_model.m          - Plots final fault-slip model
  plot_psi_posteriorPDF.m   - Plots ensemble statistics
  run_psi.sh                - Run the PSI inversion
  run_res.sh                - Run post-processing of the PSI inversion

5 INSTALATION
===================

  1) Compile Green's functions computation codes in the <dwn> folder
  2) Copy the required third-party PT and Time_2d libraries into <lib> folder
  3) Set your compilers in MAKEFILE in the <src> folder
  4) Compile the PSI by using the MAKEFILE, executable binaries should appear
  
6 EXECUTE INVERSION
===================

  1) Set input parameters into files of the <input> folder
  2) Run ./stations to prepare list of stations in X,Y coordinates (or prepare it manually)
  3) Run ./prepare to prepare input files for Green's functions computation
  4) Compute Green's functions and copy NEZsor.dat into PSI root folder
  5) Prepare observed data into NEZwave.dat file in the PSI root folder (ASCII file)
  6) Run ./waves to prepare vector of observed data var_Unez.tmp
  7) Execute PSI inversion using  ./run_psi.sh bash script
  8) Execute PSI post-processing using  ./run_res.sh bash script
  9) Plot and see results in the <inv> folder by MATLAB scripts in the PSI root folder
  
  Note: See connected example files for the structure of ASCII input files and observed data.