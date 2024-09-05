!---------------------------------------------------------------------
!
!  Initiation from INPUT files
!
!---------------------------------------------------------------------
!
!  Author: Miroslav Hallo and Frantisek Gallovic
!  Charles University, Faculty of Mathematics and Physics
!  E-mail: hallo@karel.troja.mff.cuni.cz
!  Revision 11/2017: The initial version
!
!  Copyright (C) 2017,2018  Miroslav Hallo and Frantisek Gallovic
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

MODULE InputBox
!---------------------------------------------------------------------
!  Module containing parameters from input files and others
!---------------------------------------------------------------------
    ! Constants
    REAL(8),PARAMETER:: PI=3.1415926535d0
    
    ! Direct input parameters
    INTEGER:: nfmax,NRseis,NRgps,np,NSeg,nfc
    REAL(8):: T,TS,T1,T2,artifDT,Mfix
    INTEGER,ALLOCATABLE,DIMENSION(:):: NL,NW
    REAL(8),ALLOCATABLE,DIMENSION(:):: strike,dip,rake,fc1,fc2
    REAL(8),ALLOCATABLE,DIMENSION(:):: hypodepth,leng,widt,epicW,epicL
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: offsets,lrake
    ! Prepared input parameters
    REAL(8),ALLOCATABLE,DIMENSION(:):: dL,dW,elem,stainL1,stainT
    REAL(8):: dt,df
    INTEGER:: iT1,iT2,nT,nrake
    INTEGER,ALLOCATABLE:: stainfo(:,:),fcsta(:)
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: mu,lambda
    INTEGER:: Ssvd,NSTAcomp,Nseis,Ngps,Nslip,Msvd
    
    ! Direct drive parameters
    REAL(8):: T0,thr_velo(3),thr_rise(3),thr_peak(3),thr_rake(2)
    REAL(8):: thr_hypo(3),thr_delay(3),thr_M0,thr_spline(2),ptrans
    INTEGER:: stdexist,syntdata,cdtype,MaxNSpline,debugMode,ZSCtop
    REAL(8),ALLOCATABLE,DIMENSION(:):: segDelay
    CHARACTER(len=256):: input_Uvec,input_tfin
    INTEGER,ALLOCATABLE,DIMENSION(:):: NPL,NPW
    ! Prepared drive parameters
    INTEGER:: iT0
    REAL(8),ALLOCATABLE,DIMENSION(:):: dPL,dPW
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: cellL,cellW,cpL,cpW
    ! Prepared drive ti(me)gr(id) parameters
    REAL(4):: tigrDD
    INTEGER(4),ALLOCATABLE,DIMENSION(:):: tigrNL,tigrNW
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: tigrL,tigrW
    ! PT global control
    INTEGER:: nchainsG,nstepsG,iburnG
    
    ! Direct sources.dat parameters
    REAL(4),ALLOCATABLE,DIMENSION(:,:,:):: sourNgps,sourEgps,sourZgps,strikeGPS,dipGPS
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: Rmin
!---------------------------------------------------------------------
END MODULE



SUBROUTINE InitInput()
!---------------------------------------------------------------------
!  Read parameters from input.dat and station_info.dat and prepare parameters
!---------------------------------------------------------------------
    use InputBox
    implicit none
    integer i,j,k,kk,ndepth
    real(8),allocatable,dimension(:):: depth,vp,vs,rho
    real(8):: dum,staN,staE
    
    !----------------------------------
    !  Read parameters from input file
    open(10,file='input/input.dat',action='read',status='old')
    read(10,*)
    read(10,*) nfmax
    read(10,*)
    read(10,*) np
    read(10,*)
    read(10,*) T,TS,T1,T2
    read(10,*)
    read(10,*) artifDT
    read(10,*)
    read(10,*) NRseis,NRgps
    read(10,*)
    read(10,*) Mfix
    read(10,*)
    read(10,*) NSeg
    allocate(NL(NSeg),NW(NSeg),strike(NSeg),dip(NSeg),rake(NSeg))
    allocate(hypodepth(NSeg),leng(NSeg),widt(NSeg))
    allocate(epicW(NSeg),epicL(NSeg),offsets(2,NSeg))
    read(10,*)
    read(10,*) nrake
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) (strike(i),dip(i),rake(i),i=1,NSeg)
    read(10,*)
    read(10,*) (leng(i),widt(i),i=1,NSeg)
    read(10,*)
    read(10,*) (NL(i),NW(i),i=1,NSeg)
    read(10,*)
    read(10,*) (epicL(i),epicW(i),i=1,NSeg)
    read(10,*)
    read(10,*) (hypodepth(i),i=1,NSeg)
    read(10,*)
    read(10,*) (offsets(1,i),offsets(2,i),i=1,NSeg)
    read(10,*)
    read(10,*) nfc   ! Number of frequency bands
    allocate(fc1(nfc),fc2(nfc))
    read(10,*) (fc1(i),fc2(i),i=1,nfc)
    close(10)
    
    !----------------------------------
    !  Area of spatial discretized areas
    allocate(dL(NSeg),dW(NSeg),elem(NSeg))
    dL(:)=leng(:)/dble(NL(:))
    dW(:)=widt(:)/dble(NW(:))
    elem(:)=dL(:)*dW(:)
    
    !----------------------------------
    !  Read station info from file
    ! 1-3: For N, E, Z components either 1 or 0 to select whether the station component is used
    ! 4: Refers to the number of the filter range specified in file input.dat.
    ! 5: Station name (used just for seismogram plotting)
    if(NRseis>0)then
      dt=T/dble(np)
      df=1.d0/T
      iT1=int(dble(T1)/dt+1.d0)
      iT2=int(dble(T2)/dt+1.d0)
      nT=iT2-iT1+1
      open(11,file='input/station_info.dat',action='read',status='old')
      allocate(stainfo(3,NRseis),fcsta(NRseis),stainL1(NRseis),stainT(NRseis))
      do i=1,NRseis
        read(11,*)stainfo(1:3,i),fcsta(i),stainL1(i),stainT(i)
      enddo
      close(11)
    endif
    
    !----------------------------------
    ! Read velocity model for mu and lambda
    allocate(mu(maxval(NL),maxval(NW),NSeg))
    allocate(lambda(maxval(NL),maxval(NW),NSeg))
    
    open(12,file='input/crustal.dat',action='read',status='old')
    read(12,*)
    read(12,*)
    read(12,*) ndepth
    allocate(depth(ndepth),vp(ndepth),vs(ndepth),rho(ndepth))
    read(12,*)
    read(12,*)
    do i=1,ndepth
      read(12,*) depth(i),vp(i),vs(i),rho(i)
    enddo
    close(12)
    
    !----------------------------------
    ! Evaluating mu and lambda
    do k=1,NSeg
      do i=1,NW(k)
        dum=(hypodepth(k)+(epicW(k)-dW(k)*(dble(i)-.5d0))*sin(dip(k)/180.d0*PI))/1000.d0
        if(dum>depth(ndepth))then
          mu(1:NL(k),i,k)=rho(ndepth)*vs(ndepth)**2*1.d9
          lambda(1:NL(k),i,k)=rho(ndepth)*vp(ndepth)**2*1.d9-2*mu(1:NL(k),i,k)
        else
          do j=1,ndepth
            if(dum<depth(j))exit
          enddo
          mu(1:NL(k),i,k)=rho(j-1)*vs(j-1)**2*1.d9
          lambda(1:NL(k),i,k)=rho(j-1)*vp(j-1)**2*1.d9-2*mu(1:NL(k),i,k)
        endif
      enddo
    enddo
    deallocate(depth,vp,vs,rho)
    
    !----------------------------------
    ! Process rake-angle scatter
    allocate(lrake(nrake,NSeg))
    if(nrake==1)then
      do k=1,NSeg
        lrake(1,k)=rake(k)
      enddo
    elseif(nrake==2)then
      do k=1,NSeg
        lrake(1,k)=rake(k)-45.
        if(lrake(1,k)<-180.)then
          lrake(1,k)=lrake(1,k)+360.
        endif
        lrake(2,k)=rake(k)+45.
        if(lrake(2,k)>180.)then
          lrake(2,k)=lrake(2,k)-360.
        endif
      enddo
    else
      write(*,*) ' Error! Not supported number of rakes!!'
    endif
    
    !----------------------------------
    !  Number of time samples of the slip velocity model
    if(NRseis>0)then
      Ssvd=int(TS/dt+1.d0)
      if(iT1-Ssvd<0.or.iT2>np)then
        write(*,*)' Error! Time window out of range, check input.dat!!'
        stop
      endif
    else
      dt=1.d0
      df=1.d0
      np=0
      Ssvd=1                      !Just slip value
    endif
    
    !----------------------------------
    !  Number of total used components 
    Ngps=NRgps*3
    if(NRseis>0)then
      NSTAcomp=sum(stainfo(:,:))
      Nseis=nT*NSTAcomp
    else
      NSTAcomp=0
      Nseis=0
    endif
    
    !----------------------------------
    !  Number of total synthetic data for inversion
    Msvd=sum(NW(:)*NL(:))*Ssvd*nrake
    
    
    
!---------------------------------------------------------------------
!  Read parameters from drive.dat file and prepare PSI drive parameters
!---------------------------------------------------------------------

    !----------------------------------
    !  Allocate arrays and read parameters from input file
    allocate(NPL(NSeg),NPW(NSeg),segDelay(NSeg))
    
    !  Read parameters from drive file
    open(19,file='input/drive.dat',action='read',status='old')
    read(19,*)
    read(19,*) debugMode
    read(19,*)
    read(19,*) syntdata
    read(19,*)
    read(19,*) cdtype
    read(19,*)
    read(19,*) stdexist
    read(19,*)
    read(19,*) input_Uvec
    read(19,*)
    read(19,*) input_tfin
    read(19,*)
    read(19,*)
    read(19,*) nchainsG
    read(19,*)
    read(19,*) nstepsG
    read(19,*)
    read(19,*) iburnG
    read(19,*)
    read(19,*)
    read(19,*) MaxNSpline
    read(19,*)
    read(19,*) ZSCtop
    read(19,*)
    read(19,*) tigrDD
    read(19,*)
    read(19,*) T0
    read(19,*)
    read(19,*) segDelay(1:NSeg)
    read(19,*)
    read(19,*) NPL(1:NSeg)
    read(19,*)
    read(19,*) NPW(1:NSeg)
    read(19,*)
    read(19,*)
    read(19,*) thr_velo(1:3)
    read(19,*)
    read(19,*) thr_rise(1:3)
    read(19,*)
    read(19,*) thr_peak(1:3)
    read(19,*)
    read(19,*) thr_rake(1:2)
    read(19,*)
    read(19,*) thr_hypo(1:3)
    read(19,*)
    read(19,*) thr_delay(1:3)
    read(19,*) 
    read(19,*) thr_M0
    read(19,*)
    read(19,*) thr_spline(1:2)
    read(19,*) 
    read(19,*) ptrans
    close(19)
    
    !debugMode=0   ! Run in debug mode with full log
    !syntdata=0 ! Type of the used data (0=real_data, 1=synth_data 1, 2=synth_data 2, -1=read_tfin)
    !cdtype=3 ! Type of the used covariance matrix (1=SACF, 2=ACF, 3=diagonal)
    !stdexist=1 ! Compute CD or read already computed standardized data (0=compute, 1=read)
    !input_Uvec='var_Unez.tmp' ! if syntdata=0 (real data)
    !input_tfin='var_tfin.tmp' ! if syntdata=-1 (tfin)
    !-------------
    !nchainsG = 30    ! Define number of chains
    !nstepsG = 100000 ! Number of chain steps per temperature
    !iburnG = 500      ! Number of burn in samples for skipping models
    !-------------
    !MaxNSpline=30 ! Maximal number of efficient spline control points
    !ZSCtop=0 ! Slip constraint on top of faults (0=zero slip constraint on top, 1=no constraint on top)
    !tigrDD = 500. ! Spacing (meters) among points for travel-time computation
    !T0=-0.8d0 ! Time shift of the synthetic seismograms [sec] (negative means later synthetics)
    !segDelay=(/ 0.8d0 /) ! Segment time delay [sec]
    !NPL=(/ 9 /) ! Number of control points along strike
    !NPW=(/ 7 /) ! Number of control points along dip
    !-------------
    !thr_velo=(/ 1500.d0, 5000.d0, 50.d0 /) ! Min, Max, stepSigma
    !thr_rise=(/ 1.0d0, 7.0d0, 0.05d0 /) ! Min, Max, stepSigma
    !thr_peak=(/ 0.4d0, 2.8d0, 0.05d0 /) ! Min, Max, stepSigma
    !thr_rake=(/ 45.d0, 1.d0 /) ! +/-Deg_max, stepSigma
    !thr_hypo=(/ 1000.d0, 2000.d0, 50.d0/) ! +/-L_max, +/-W_max, stepSigma
    !thr_delay=(/ 0.8d0, 0.0d0, 0.05d0/) !  +/-delay_max_delta (1: 1st segment, 2: other segments, 3: stepSigma) (segment delay has to be >= 0)
    !thr_M0= 2.d18 ! +/-M0_max
    !thr_spline=(/ 50.d0, 0.01d0 /) ! pos-stepSigma, val-stepSigma
    !ptrans=0.05 ! Probability to change dimension (*100%)
    
    
    !----------------------------------
    !  Additional temporal shift between synth a reald data
    iT0=int(T0/dt)
    !  Check for error
    if(NRseis>0)then
      Ssvd=int(TS/dt+1.d0)       !number of time samples of the slip velocity model
      if(iT1-iT0-Ssvd<0.or.iT2-iT0>np)then
        write(*,*)' Error! Time window out of range, check input...'
        stop
      endif
    else
      dt=1.d0
      df=1.d0
      np=0
      Ssvd=1                      !Just slip value
    endif
    
    !----------------------------------
    !  Generates cell grid
    allocate(cellL(maxval(NL),NSeg),cellW(maxval(NW),NSeg))
    do kk=1,NSeg
      do i=1,NL(kk)
        cellL(i,kk)=(dble(i)-.5d0)*dL(kk)
      enddo
      do j=1,NW(kk)
        cellW(j,kk)=(dble(j)-.5d0)*dW(kk)
      enddo
    enddo
    
    !----------------------------------
    !  Spacing between control points
    allocate(dPL(NSeg),dPW(NSeg))
    dPL(:)=leng(:)/dble(NPL(:)-1)
    dPW(:)=widt(:)/dble(NPW(:)-1)
    
    !----------------------------------
    !  Generates control points grid
    allocate(cpL(maxval(NPL),NSeg),cpW(maxval(NPW),NSeg))
    do kk=1,NSeg
      do i=1,NPL(kk)
        cpL(i,kk)=(dble(i-1))*dPL(kk)
      enddo
      do j=1,NPW(kk)
        cpW(j,kk)=(dble(j-1))*dPW(kk)
      enddo
    enddo
    
    !----------------------------------
    !  Generates dense grid for rupture time computation
    !        ti(me)gr(id) parameters
    allocate(tigrNL(NSeg),tigrNW(NSeg))
    do kk=1,NSeg
      tigrNL(kk)=int(ceiling(real(leng(kk))/tigrDD))+1
      tigrNW(kk)=int(ceiling(real(widt(kk))/tigrDD))+1
    enddo
    
    allocate(tigrL(maxval(tigrNL),NSeg),tigrW(maxval(tigrNW),NSeg))
    do kk=1,NSeg
      do i=1,tigrNL(kk)
        tigrL(i,kk)=(dble(i-1))*dble(tigrDD)
      enddo
      do j=1,tigrNW(kk)
        tigrW(j,kk)=(dble(j-1))*dble(tigrDD)
      enddo
    enddo
    
    
    
!---------------------------------------------------------------------
!  Read parameters from sources.dat and stations.dat files
!     (needed for distance-dependent weights -> CD matrix calculation)
!---------------------------------------------------------------------

    !----------------------------------
    !  Allocate
    allocate(sourNgps(maxval(NL),maxval(NW),NSeg),sourEgps(maxval(NL),maxval(NW),NSeg),sourZgps(maxval(NL),maxval(NW),NSeg))
    allocate(strikeGPS(maxval(NL),maxval(NW),NSeg),dipGPS(maxval(NL),maxval(NW),NSeg))
    allocate(Rmin(NRseis,NSeg))
    
    !----------------------------------
    !  Read parameters from sources.dat file
    open(13,file='sources.dat',action='read',status='old')
    do k=1,NSeg
      do j=1,NW(k)
        do i=1,NL(k)
          read(13,*)dum,sourNgps(i,j,k),sourEgps(i,j,k),sourZgps(i,j,k),strikeGPS(i,j,k),dipGPS(i,j,k)
        enddo
      enddo
    enddo
    close(13)
    sourNgps=sourNgps*1.e3
    sourEgps=sourEgps*1.e3
    sourZgps=sourZgps*1.e3

    !----------------------------------
    !  Read parameters from stations.dat file 
    !       and compute minimal distances from the leading segment
    open(14,file='stations.dat',action='read',status='old')
    do i=1,NRseis
      read(14,*) staN,staE
      do k=1,NSeg
        Rmin(i,k)=minval(sqrt((sourNgps(1:NL(k),1:NW(k),k)-staN*1.e3)**2+(sourEgps(1:NL(k),1:NW(k),k)-staE*1.e3)**2))
      enddo
    enddo
    close(14)
    
    
    
!---------------------------------------------------------------------
END SUBROUTINE


