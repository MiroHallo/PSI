!---------------------------------------------------------------------
!
!  Parametric Slip Inversion (PSI)
!
!  Hallo, M., Gallovic, F. (2020): Bayesian self-adapting fault slip inversion 
!   with Green's functions uncertainty and application on the 2016 Mw7.1 Kumamoto
!   earthquake, Journal of Geophysical Research: Solid Earth, 125, e2019JB018703.
!
!---------------------------------------------------------------------
!
!  Author: Miroslav Hallo
!  Charles University, Faculty of Mathematics and Physics
!  E-mail: hallo@karel.troja.mff.cuni.cz
!  Revision 11/2017: The initial version
!
!  Copyright (C) 2017  Miroslav Hallo
!
!  This program is published under the GNU General Public License (GNU GPL).
!  This program is free software: you can modify it and/or redistribute it
!  or any derivative version under the terms of the GNU General Public
!  License as published by the Free Software Foundation, either version 3
!  of the License, or (at your option) any later version.
!  This code is distributed in the hope that it will be useful, but WITHOUT
!  ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
!  and don't remove their names from the code.
!  You should have received copy of the GNU General Public License along
!  with this program. If not, see <http://www.gnu.org/licenses/>.
!
!---------------------------------------------------------------------

PROGRAM  psi
!---------------------------------------------------------------------
!  Main program of the Parametric Slip Inversion (PSI), Parallel Tempering
!---------------------------------------------------------------------
    implicit none
    
    !----------------------------------
    write(*,*) ""
    write(*,*) "------------------------------------------------------"
    write(*,*) " ._       _.      _____   _____ _____"
    write(*,*) " |\ \\_// /|     |  __ \ / ____|_   _|"
    write(*,*) " \/       \/     | |__) | (___   | |"
    write(*,*) "_|_ V   V _|_    |  ___/ \___ \  | |"
    write(*,*) "= - = * = - =    | |     ____) |_| |_"
    write(*,*) "   \_____/       |_|    |_____/|_____|"
    write(*,*) "    /_ _\      Parametric Slip Inversion"
    write(*,*) "   |- W -\--..      - Parallel Tempering"
    write(*,*) "  | |   | |   \"
    write(*,*) "   \|_ _|/ \.  \    Miroslav Hallo"
    write(*,*) "  /__\ /__\  \__>   2017-2019"
    
    !----------------------------------
    ! Read input files and init input variables
    write(*,*) "------------------------------------------------------"
    write(*,*) "Initiation of the input"
    call InitInput()
    
    !----------------------------------
    ! Initiate core matrices H and Uvec
    write(*,*) "------------------------------------------------------"
    write(*,*) "Initiation of the core"
    call InitCore()
    
    !----------------------------------
    ! Execute of the inversion
    write(*,*) "------------------------------------------------------"
    write(*,*) "Execute of the inversion"
    call InvExe()
    
    write(*,*) "PSI DONE"
!---------------------------------------------------------------------
END PROGRAM



    
    
    
    