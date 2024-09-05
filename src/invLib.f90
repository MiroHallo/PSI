!---------------------------------------------------------------------
!
!  Module with functions used by PSI
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

MODULE LibBox
!---------------------------------------------------------------------
!  Module containing global functions
!---------------------------------------------------------------------
CONTAINS



FUNCTION RYoffe(N,dt,t0,tR,tP,D)
!---------------------------------------------------------------------
!  Generates discrete time array with Regularized Yoffe function
!  (Result in RYoffe is slip-rate time function)
!  AUTHOR: Miroslav Hallo, 2017
!  Following: Tinti E., Fukuyama E., Piatanesi A., Cocco M. (2005). A Kinematic Source-Time
!             Function Compatible with Earthquake Dynamics, BSSA, 95(4), 1211–1223.
!  INPUT: N .. number of samples
!         dt .. sampling interval [sec]
!         t0 .. rupture time [sec] when the 1st sample is 0.0 time
!         tR .. rise time [sec] relatively to t0
!         tP .. peak time [sec] relatively to t0
!         D .. total slip [m]
!---------------------------------------------------------------------
    implicit none
    real(8),parameter:: pi=3.1415926535d0
    integer,intent(in):: N
    real(8),intent(in):: dt,t0,tR,tP,D
    real(8):: C1,C2,C3,C4,C5,C6,t,tS,K,tfc(N),RYoffe(N)
    integer:: i
    !----------------------------------
    !  Prepare tS
    tS=tP/1.27d0
    if(tR<2*tS)then
      !write(*,*)' Error! Risetime < 2*Smootime (Yoffe not defined)!!!'
      !stop
      RYoffe(1)=-12345.d0
      RETURN
    endif
    
    !----------------------------------
    !  Equation A13
    K=D*2.d0/(pi*tR*tS**2)
    RYoffe=0.d0
    tfc=dt*(/(dble(i),i=0,N-1)/)-t0
    do i=1,N
      t=tfc(i)
      if(t<0)then
          cycle
      elseif(t<tS)then
          C1=(.5d0*t+.25d0*tR)*sqrt(t*(tR-t))+(t*tR-tR**2)* &
             dasin(sqrt(t/tR))-(3.d0/4.d0*tR**2)*datan(sqrt((tR-t)/t))
          C2=(3.d0/8.d0)*pi*tR**2
          RYoffe(i)=K*(C1+C2)
      elseif(t<2.d0*tS)then
          C1=(.5d0*t+.25d0*tR)*sqrt(t*(tR-t))+(t*tR-tR**2)* &
             dasin(sqrt(t/tR))-(3.d0/4.d0*tR**2)*datan(sqrt((tR-t)/t))
          C2=(3.d0/8.d0)*pi*tR**2
          C3=(tS-t-.5d0*tR)*sqrt((t-tS)*(tR-t+tS))+ &
             tR*2.d0*(tR-t+tS)*dasin(sqrt((t-tS)/tR))+ &
             (3.d0/2.d0)*(tR**2)*datan(sqrt((tR-t+tS)/(t-tS)))
          RYoffe(i)=K*(C1-C2+C3)
      elseif(t<tR)then
          C1=(.5d0*t+.25d0*tR)*sqrt(t*(tR-t))+(t*tR-tR**2)* &
             dasin(sqrt(t/tR))-(3.d0/4.d0*tR**2)*datan(sqrt((tR-t)/t))
          C3=(tS-t-.5d0*tR)*sqrt((t-tS)*(tR-t+tS))+ &
             tR*2.d0*(tR-t+tS)*dasin(sqrt((t-tS)/tR))+ &
             (3.d0/2.d0)*(tR**2)*datan(sqrt((tR-t+tS)/(t-tS)))
          C4=(-tS+.5d0*t+.25d0*tR)*sqrt((t-2.d0*tS)*(tR-t+2.d0*tS)) - &
             tR*(tR-t+2.d0*tS)*dasin(sqrt((t-2.d0*tS)/tR)) - &
            (3.d0/4.d0)*(tR**2)*datan(sqrt((tR-t+2.d0*tS)/(t-2.d0*tS)))
          RYoffe(i)=K*(C1+C3+C4)
      elseif(t<tR+tS)then
          C3=(tS-t-.5d0*tR)*sqrt((t-tS)*(tR-t+tS))+ &
             tR*2.d0*(tR-t+tS)*dasin(sqrt((t-tS)/tR))+ &
             (3.d0/2.d0)*(tR**2)*datan(sqrt((tR-t+tS)/(t-tS)))
          C4=(-tS+.5d0*t+.25d0*tR)*sqrt((t-2.d0*tS)*(tR-t+2.d0*tS)) - &
             tR*(tR-t+2.d0*tS)*dasin(sqrt((t-2.d0*tS)/tR)) - &
            (3.d0/4.d0)*(tR**2)*datan(sqrt((tR-t+2.d0*tS)/(t-2.d0*tS)))
          C5=(pi/2.d0)*tR*(t-tR)
          RYoffe(i)=K*(C5+C3+C4)
      elseif(t<tR+2.d0*tS)then
          C4=(-tS+.5d0*t+.25d0*tR)*sqrt((t-2.d0*tS)*(tR-t+2.d0*tS)) - &
             tR*(tR-t+2.d0*tS)*dasin(sqrt((t-2.d0*tS)/tR)) - &
            (3.d0/4.d0)*(tR**2)*datan(sqrt((tR-t+2.d0*tS)/(t-2.d0*tS)))
          C6=(pi/2.d0)*tR*(2.d0*tS-t+tR)
          RYoffe(i)=K*(C4+C6)
      else
          exit
      endif
    enddo
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION binarysearch(length,array,value)
!---------------------------------------------------------------------
! Given an array and a value, returns the index of the element that
! is closest to, but less than, the given value.
! Uses a binary search algorithm.
! Used by interp2D function
!---------------------------------------------------------------------
    implicit none
    real(8),parameter:: dpr=1.d-9
    integer,intent(in):: length
    real(8),dimension(length),intent(in):: array
    real(8),intent(in):: value
    integer:: binarysearch
    integer:: left, middle, right
    !----------------------------------
    !  Search for index
    left = 1
    right = length
    do
      if (left > right) then
        exit
      endif
      middle = nint((left+right) / 2.0)
      if ( abs(array(middle) - value) <= dpr) then
        binarySearch = middle
        RETURN
      elseif (array(middle) > value) then
        right = middle - 1
      else
        left = middle + 1
      endif
    enddo
    binarysearch = right
    RETURN
!---------------------------------------------------------------------
END FUNCTION



REAL(8) FUNCTION interp2D(x_len,x_array,y_len,y_array,f,x,y)
!---------------------------------------------------------------------
! This function uses bilinear interpolation to estimate the value
! of a function f at point (x,y)
! f is assumed to be sampled on a regular grid, with the grid x values specified
! by x_array and the grid y values specified by y_array
! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: x_len,y_len           
    real(8),dimension(x_len),intent(in):: x_array
    real(8),dimension(y_len),intent(in):: y_array
    real(8),dimension(x_len,y_len),intent(in):: f
    real(8),intent(in):: x,y
    real(8):: denom,x1,x2,y1,y2
    integer:: i,j
    !----------------------------------
    !  Search for index
    i = binarysearch(x_len, x_array, x)
    j = binarysearch(y_len, y_array, y)
    !----------------------------------
    !  Search for index
    x1 = x_array(i)
    x2 = x_array(i+1)
    y1 = y_array(j)
    y2 = y_array(j+1)
    denom = (x2 - x1)*(y2 - y1)
    interp2D = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
       f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION bspline(Np,Nx,Ny,X,Y,Z,XI,YI)
!---------------------------------------------------------------------
! 2D Biharmonic spline interpolation
!   Sandwell, D. T. (1987), Biharmonic spline interpolation of GEOS-3
!   ans SEASAT altimeter data, Geophysical Research Letters, Vol. 2, 
!   p. 139 – 142.
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Np,Nx,Ny
    real(8),intent(in):: X(Np),Y(Np),Z(Np),XI(Nx),YI(Ny)
    real(8):: bspline(Nx,Ny)
    real(8):: GG(Np,Np),g(Np),magx,m(Np)
    integer:: i,j,k,IPIV(Np),INFO
    !----------------------------------
    ! Compute GG matrix for GG*m = Z inversion problem
    GG=0.d0
    do i=1,Np
      do j=1,i-1
        magx=sqrt((X(i)-X(j))**2 + (Y(i)-Y(j))**2)
        if(magx>=1.d-7)then
          GG(i,j)=(magx**2)*(log(magx)-1.d0)
          GG(j,i)=GG(i,j)
        endif
      enddo
    enddo
    !----------------------------------
    ! Compute model "m" where data "GG" is equal to "Z" (m = GG\Z)
    ! LAPACK linear solver (use -llapack -lblas for gfortran)
    ! LAPACK linear solver (use -mkl -DMKL for ifort)
    m=Z
    call DGESV(Np,1,GG,Np,IPIV,m,Np,INFO)
    if(INFO.ne.0)then
      bspline(1,1)=-12345.d0
      RETURN
    endif
    !----------------------------------
    ! Find 2D interpolated surface ZI through XI, YI grid points
    do i=1,Nx
      do j=1,Ny
        do k=1,Np
          magx=sqrt((XI(i)-X(k))**2 + (YI(j)-Y(k))**2)
          if(magx>=1.d-7)then
            g(k)=(magx**2)*(log(magx)-1.d0)
          else
            g(k)=(magx**2)*(-100.d0)
          endif
        enddo
        bspline(i,j) = sum(g*m)
      enddo
    enddo
    RETURN
!---------------------------------------------------------------------
END FUNCTION



!---------------------------------------------------------------------
END MODULE


