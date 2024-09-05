!---------------------------------------------------------------------
!
!  Initiation core matrices and variables
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

MODULE CoreBox
!---------------------------------------------------------------------
!  Module containing global core variables
!---------------------------------------------------------------------
#ifdef REAL4
    REAL(4),ALLOCATABLE:: H(:,:),Uvec(:),tfin(:),CD4(:,:)
#else
    REAL(8),ALLOCATABLE:: H(:,:),Uvec(:),tfin(:)
#endif
    REAL(8),ALLOCATABLE:: CD(:,:)
!---------------------------------------------------------------------
END MODULE



SUBROUTINE InitCore()
!---------------------------------------------------------------------
!  Creates matrices G of elementary synthetics (from dwn), data vector
!  and other global variables which can be prepared prior inversion
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    implicit none
#ifdef REAL4
    real(4),allocatable:: Hnew(:,:),Uvec_tmp(:)
#else
    real(8),allocatable:: Hnew(:,:),Uvec_tmp(:)
#endif
    real(8) dum
    integer i,j
    
    !----------------------------------
    ! Print some input variables
    if(debugMode>0)then
      write(*,*) " Area of spatial discretized areas [km^2] all segments:"
      write(*,'(2F10.4)') elem/1000000.d0
      write(*,*) " Spatial discretized size along strike [km] all segments:"
      write(*,'(2F10.4)') dL/1000.d0
      write(*,*) " Spatial discretized size along dip [km] all segments:"
      write(*,'(2F10.4)') dW/1000.d0
      write(*,*) " Number of time samples of one waveform:"
      write(*,'(I6)') nT
      write(*,*) " Whole data vector length [samples]:"
      write(*,'(I8)') Nseis
      write(*,*) " Source-time function length [samples]:"
      write(*,'(I6)') Ssvd
      write(*,*) " Synth. dt [sec]:"
      write(*,'(F10.4)') dt
      write(*,*) " Synth. df [Hz]:"
      write(*,'(F10.4)') df
      write(*,*) " Number of receivers (waveforms, GPS):"
      write(*,'(2I5)') NRseis,NRgps
      write(*,*) " Component usage (N E Z)  Filter_freq [Hz]"
      do i=1,NRseis
        write(*,'(3I3,2F8.4)') stainfo(:,i),fc1(fcsta(i)),fc2(fcsta(i))
      enddo
      write(*,*) " First computed Lame mu parameter"
      write(*,'(E15.5)') mu(1,1,1)
      write(*,*) " Number of computed rakes"
      write(*,'(I5)') nrake
      write(*,*) " Used rakes"
      do i=1,NSeg
        write(*,'(2F9.1)') lrake(:,i)
      enddo
      ! Print some drive variables
      write(*,*) " Type of data used (0..real(file), -1..tfin(file), >1..synt.gen:"
      write(*,'(I5)') syntdata
      write(*,*) " Additional temporal shift [samples]:"
      write(*,'(I5)') iT0
    endif
    
    !----------------------------------
    ! Allocate
    allocate(H(Nseis,Msvd),tfin(Msvd),Uvec(Nseis))
    
    !----------------------------------
    ! Init GPS
    if(NRgps>0)then
      write(*,*)' GPS data are not supported at this version of program'
    endif
    
    !----------------------------------
    ! Read standardized data from files if possible (already computed by first run)
    if(stdexist>0)then
      write(*,*)' Read standardized data from var_Unezstd.tmp'
      open(421,file='var_Unezstd.tmp')
      do i=1,nT
        read(421,*) dum,(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(421)
      
      write(*,*)' Read standardized data from var_Hstd.tmp'
      open(422,FILE='var_Hstd.tmp',ACCESS='stream',FORM='unformatted')
      do i=1,Nseis
        read(422) H(i,1:Msvd)
      enddo
      close(422)
      
    !----------------------------------
    ! First run (read data, compute CD and create standardized data)
    else
      ! Creates matrices H of elementary synthetics in time domain
      write(*,*)' Creating matrix H'
      call CreateH()
      
      ! Write full H to file
      open(233,FILE='var_H.tmp',ACCESS='stream',FORM='unformatted')
      do i=1,Nseis
        write(233) H(i,1:Msvd)
      enddo
      close(233)
      
      ! Write for debugging
      if(debugMode>1)then
        open(231,form='formatted',file='tempH.tmp')
          do i=1,Nseis
            write(231,'(5E13.5)') H(i,1:5)
          enddo
        close(231)
      endif
      
      ! Creates vector Uvec of data
      write(*,*)' Creating data vector Uvec'
      call CreateUvec()
      
      ! Prepare covariance matrix following Hallo and Gallovic (2016)
      write(*,*)' Creating covariance matrix CD'
      allocate(CD(Nseis,Nseis))
#ifdef REAL4
      allocate(CD4(Nseis,Nseis))
#endif
      call CreateCD()
      
      ! Prepare standardized Uvec
      write(*,*)' Creating standardized data from Uvec'
#ifdef MKL
      allocate(Uvec_tmp(Nseis))
      Uvec_tmp=Uvec
#ifdef REAL4
      call sgemv('N',Nseis,Nseis,1.,CD4,Nseis,Uvec_tmp,1,0.,Uvec,1)
#else
      call dgemv('N',Nseis,Nseis,1.d0,CD,Nseis,Uvec_tmp,1,0.d0,Uvec,1)
#endif
      deallocate(Uvec_tmp)
#else
#ifdef REAL4
      Uvec=matmul(CD4,Uvec)
#else
      Uvec=matmul(CD,Uvec)
#endif
#endif
    
      ! Write standardized waveforms
      open(232,FILE='var_Unezstd.tmp')
      do i=1,nT
        write(232,'(1000E17.7E3)')dt*(iT1-1+i-1),(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(232)
    
      ! Prepare standardized H
      write(*,*)' Creating standardized data from H'
      allocate(Hnew(Nseis,Msvd))
#ifdef MKL
#ifdef REAL4
      call sgemm('N','N',Nseis,Msvd,Nseis,1.,CD4,Nseis,H,Nseis,0.,Hnew,Nseis)
#else
      call dgemm('N','N',Nseis,Msvd,Nseis,1.d0,CD,Nseis,H,Nseis,0.d0,Hnew,Nseis)
#endif
#else
#ifdef REAL4      
      Hnew=matmul(CD4,H)
#else
      Hnew=matmul(CD,H)
#endif
#endif
      H=Hnew
      deallocate(Hnew)

      ! Write full H to file
      open(233,FILE='var_Hstd.tmp',ACCESS='stream',FORM='unformatted')
      do i=1,Nseis
        write(233) H(i,1:Msvd)
      enddo
      close(233)
      
      ! Write H for debugging
      if(debugMode>1)then
        open(231,form='formatted',file='tempHstd.tmp')
        do i=1,Nseis
          write(231,'(5E13.5)') H(i,1:5)
        enddo
        close(231)
      endif
    
      ! Deallocate (we dont need it any more)
      deallocate(CD)
#ifdef REAL4
      deallocate(CD4)
#endif
      
    endif
    
!---------------------------------------------------------------------
CONTAINS



SUBROUTINE CreateH()
!---------------------------------------------------------------------
!  Creates matrices of elementary synthetics (from dwn)
!     IFFT of dwn into time domain, filtration, integration, time shifts
!---------------------------------------------------------------------
    implicit none
    real*4 dumarr(6),dtr4,f1r4,f2r4
    integer dumi
    complex*16,dimension(:,:,:),allocatable:: cirN,cirE,cirZ
    complex*16,dimension(:),allocatable:: cseis
    real*4,dimension(:),allocatable:: fltr4
    integer i,j,jj,k,kk,r,mm,jw,jl,SegShift,RakShift
    
    open(20,form='unformatted',file='NEZsor.dat')
    !if problems with NEZsor.dat file appears, check that resort.f90 writes unformatted file (and not binary)!
    if(debugMode>0)then
      write(*,'(A,F9.3,A4)')'  Correct GFs for artifical time delay ',artifDT,' sec.'
    endif
    
    RakShift=sum(NW(:)*NL(:))*Ssvd
    do r=1,nrake
      do kk=1,NSeg
        allocate(cirN(min(nfmax,np),NL(kk)*NW(kk),NRseis))
        allocate(cirE(min(nfmax,np),NL(kk)*NW(kk),NRseis))
        allocate(cirZ(min(nfmax,np),NL(kk)*NW(kk),NRseis))
        SegShift=((r-1)*RakShift)+sum(NW(1:kk-1)*NL(1:kk-1))*Ssvd
        do i=1,NRseis
          j=0
          read(20) dumi
          do jw=1,NW(kk)
            do jl=1,NL(kk)
              j=j+1
              read(20) dumi
              do k=1,nfmax
                read(20) dumarr
                if(k>np)cycle
                cirN(k,j,i)=cmplx(dumarr(1),dumarr(4))*mu(jl,jw,kk)*elem(kk)*dt
                cirE(k,j,i)=cmplx(dumarr(2),dumarr(5))*mu(jl,jw,kk)*elem(kk)*dt
                cirZ(k,j,i)=cmplx(dumarr(3),dumarr(6))*mu(jl,jw,kk)*elem(kk)*dt
              enddo
            enddo
          enddo
        enddo
      
        dtr4=real(dt)
        j=0
        do jj=1,NRseis
          f1r4=real(fc1(fcsta(jj)));f2r4=real(fc2(fcsta(jj)))
          do k=1,3
            if(stainfo(k,jj)==0)cycle
            j=j+1
!$OMP parallel private(i,mm,cseis,fltr4) DEFAULT(SHARED)
            allocate(cseis(np),fltr4(np))
!$OMP do SCHEDULE(DYNAMIC,1)
            do i=1,NL(kk)*NW(kk)
              cseis=0.d0
              select case (k)
              case(1)
                cseis(1:min(nfmax,np))=cirN(1:min(nfmax,np),i,jj)*df
              case(2)
                cseis(1:min(nfmax,np))=cirE(1:min(nfmax,np),i,jj)*df
              case(3)
                cseis(1:min(nfmax,np))=cirZ(1:min(nfmax,np),i,jj)*df
              end select
              cseis(np/2+2:np)=conjg(cseis(np/2:2:-1))
              call four1(cseis,np,1)
              fltr4=real(cseis*dt)   !dt for time integration
              do mm=1,int(artifDT/dt)  ! REMOVING DWN ARTIFACT FROM THE SEISMOGRAM BEGINNING
                fltr4(mm)=fltr4(mm)*real(cos((dt*dble(mm-1)-artifDT)/artifDT*PI)/2.d0+0.5d0)
              enddo
              if(f1r4>0.)then
                CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'BP', f1r4, f2r4, dtr4, 1, np)
              else
                CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'LP', f1r4, f2r4, dtr4, 1, np)
              endif
              do mm=2,np   !time integration (velocity to displacement)
                fltr4(mm)=fltr4(mm)+fltr4(mm-1)
              enddo
              ! Write to file for debugging
              if(debugMode>1)then
                if(kk==1.and.i==1.and.jj==1.and.k==1)then
                  open(230,file='tempSynt.tmp')
                  do mm=1,np
                    write(230,*)dt*(mm-1),fltr4(mm)
                  enddo
                  close(230)
                endif
              endif
              do mm=1,nT
#ifdef REAL4
                  H((j-1)*nT+mm,SegShift+(i-1)*Ssvd+1:SegShift+i*Ssvd)=fltr4(iT1-1+mm-iT0:iT1+mm-Ssvd-iT0:-1)
#else
                  H((j-1)*nT+mm,SegShift+(i-1)*Ssvd+1:SegShift+i*Ssvd)=dble(fltr4(iT1-1+mm-iT0:iT1+mm-Ssvd-iT0:-1))
#endif
              enddo
            enddo
!$OMP end do
            deallocate(cseis,fltr4)
!$OMP end parallel
          enddo
        enddo
        deallocate(cirN,cirE,cirZ)
      enddo
    enddo
    close(20)
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE CreateUvec()
!---------------------------------------------------------------------
!  Creates vector of data (from real data or from synthetics)
!---------------------------------------------------------------------
    implicit none
    integer i,j
    real(8) dum
    
    if(syntdata==0)then
      ! Read Uvec from real data
      if(debugMode>0)then
        write(*,*)' Read Uvec from ',trim(input_Uvec)
      endif
      open(13,file=trim(input_Uvec))
      do i=1,nT
        read(13,*) dum,(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(13)
      
    elseif(syntdata>0)then
      if(debugMode>0)then
        write(*,*)' Generate synthetic slip-rate model'
      endif
      ! Create synthetic slip-rate time function --> tfin
      if(syntdata==1)then
        call synth_tfin1()
      elseif(syntdata==2)then
        call synth_tfin2()
      endif
      ! Create Uvec from STF tfin
#ifdef MKL
#ifdef REAL4
      call sgemv('N',Nseis,Msvd,1.,H,Nseis,tfin,1,0.,Uvec,1)
#else
      call dgemv('N',Nseis,Msvd,1.d0,H,Nseis,tfin,1,0.d0,Uvec,1)
#endif
#else
      Uvec=matmul(H,tfin)
#endif
      ! Write waveforms
      open(13,FILE='var_Unez.tmp')
      do i=1,nT
        write(13,'(1000E17.7E3)')dt*(iT1-1+i-1),(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(13)
      
    elseif(syntdata<0)then
      if(debugMode>0)then
        write(*,*)' Read custom target model from ',trim(input_tfin)
      endif
      open(12,file=trim(input_tfin))
      do i=1,Msvd
        read(12,*)tfin(i)
      enddo
      close(12)
      ! Create Uvec from STF tfin
#ifdef MKL
#ifdef REAL4
      call sgemv('N',Nseis,Msvd,1.,H,Nseis,tfin,1,0.,Uvec,1)
#else
      call dgemv('N',Nseis,Msvd,1.d0,H,Nseis,tfin,1,0.d0,Uvec,1)
#endif
#else
      Uvec=matmul(H,tfin)
#endif
      ! Write waveforms
      open(13,FILE='var_Unez.tmp')
      do i=1,nT
        write(13,'(1000E17.7E3)')dt*(iT1-1+i-1),(Uvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(13)
      
    else
      write(*,*)' Error! Not known type of data!!!'
      stop
    endif
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE CreateCD()
!---------------------------------------------------------------------
!  Creates covariance matrix CD from data vector Uvec 
!     following approach by Hallo and Gallovic, 2016
!---------------------------------------------------------------------
    implicit none
    real(8):: fvec(nT),CDsub(nT,nT),maxdiagCD
    real(8):: L1_low,L1_fac,d_err,wl
    real(8),allocatable:: CDinv(:,:)
    integer,allocatable,dimension(:):: L1_all
    integer:: L1,Tef
    integer:: i,j,j0,k,jj,erro
    
    ! Hard-code parameters
    L1_low = 0.8d0 ! approx. 1/3 of the minimum wave period
    L1_fac = 10000.d0 ! 15% velocity model unc.
    d_err = 0.01d0 ! Sigma to diagonal covariance matrix (cdtype=3)
    wl = 25.d0 ! Percentage of the maximum value in CD as water level on diagonal
    
    ! Allocate and prepare
    allocate(L1_all(NRseis))
    CD=0.d0
    j=0
    j0=0
    
    ! Do for all receivers
    do jj=1,NRseis
    
      ! Prepare parameter L1
      if(stainL1(jj).le.0.d0)then
        L1 = int( max(2.d0,0.5d0+max(L1_low,minval(Rmin(jj,1:NSeg))/L1_fac)/dt) )
      else
        L1 = int( max(2.d0,0.5d0+stainL1(jj)/dt) )
      endif
      L1_all(jj)=L1
      ! Prepare parameter T
      if(stainT(jj).le.0.d0)then
        Tef = int(nT/2.d0)
      else
        Tef = int( max(1.d0,0.5d0+stainT(jj)/dt) )
      endif
      
      ! Create CD matrix blocks for this receiver
      maxdiagCD=0.d0
      do k=1,3
        if(stainfo(k,jj)==0)cycle
        j = j+1
        fvec(1:nT) = Uvec((j-1)*nT+1:j*nT)
        ! Compute covariance matrix
        if(cdtype.eq.1)then
          ! SACF matrix
          call saxcf(nT,fvec,fvec,L1,0,Tef,CDsub)
        elseif(cdtype.eq.2)then
          ! ACF matrix
          call axcf(nT,fvec,fvec,L1,0,CDsub)
        else
          ! Diagonal matrix
          CDsub=0.d0
          do i=1,nT
            CDsub(i,i) = d_err**2
          enddo
        endif
        ! Compose into block-structure CD matrix
        CD((j-1)*nT+1:j*nT,(j-1)*nT+1:j*nT) = CDsub(1:nT,1:nT)
        ! Maximum of diagonal
        do i=1,nT
          if(CDsub(i,i)>maxdiagCD)maxdiagCD=CDsub(i,i)
        enddo
      enddo
      
      ! Add to diagonal station water level
      if(cdtype.ne.3)then
        do k=1,3
          if(stainfo(k,jj)==0)cycle
          j0 = j0+1
          do i=(j0-1)*nT+1,j0*nT
            CD(i,i) = CD(i,i) + wl*(maxdiagCD/100.d0)
          enddo
        enddo
      endif
    enddo
    
    ! Write L1 for debugging
    if(debugMode>0)then
      write(*,*) " L1 for used receivers in samples:"
      write(*,'(1000I6)') L1_all(1:NRseis)
      write(*,*) " L1 for used receivers in sec:"
      write(*,'(1000F6.2)') L1_all(1:NRseis)*dt
    endif
    
    ! Write to file for debugging
    if(debugMode>1)then
      ! Save L1 parameters and Cov matrix (binary)
      open(357,FILE='tempCD.tmp',ACCESS='stream',FORM='unformatted')
      write(357) int(NRseis,4)
      write(357) int(nT,4)
      write(357) int(Nseis,4)
      write(357) int(L1_all(1:NRseis),4)
      write(357) real(dble(L1_all(1:NRseis))*dt,4)
      do i=1,Nseis
        write(357) real(CD(i,1:Nseis),4)
      enddo
      close(357)
      
      ! Save L1 parameters and Cov matrix (ASCII)
      !open(357,FILE='tempCD.tmp')
      !write(357,'(3I7)') NRseis,nT,Nseis
      !write(357,'(100000I7)') L1_all(1:NRseis)
      !write(357,'(100000F7.2)') dble(L1_all(1:NRseis))*dt
      !do i=1,Nseis
        !write(357,'(100000E16.7E3)') CD(i,1:Nseis)
      !enddo
      !close(357)
    endif
    
    ! Inverse of Cholesky decomposition of CD
    do i=1,Nseis-1 
      CD(i+1:Nseis,i)=CD(i,i+1:Nseis)     ! Copy upper triangle of CD to lower
    enddo
#ifdef MKL
    call dpotrf('U',Nseis,CD,Nseis,erro)  ! Upper triangle of CD becomes U
    if(erro>0)then
      write(*,*) ' Error! The matrix CD is not positive definite!'
      stop
    endif
    call dpotri('U',Nseis,CD,Nseis,erro)  ! Inverse of CD
    call dpotrf('U',Nseis,CD,Nseis,erro)  ! Upper triangle of CD becomes U
    do i=1,Nseis-1                        
      CD(i+1:Nseis,i)=0.d0                ! Clear lower triangle
    enddo
#else
    write(*,*) ' Warning! Inverse of Cholesky decomposition by NumRec'
    write(*,*) '          It is slow and less precise! Use MKL library if you can!'
    allocate(CDinv(Nseis,Nseis))
    CDinv=CD
    call choldcsl(Nseis,CDinv,CD)
    do i=1,Nseis-1
      CD(i,i+1:Nseis)=CD(i+1:Nseis,i)     ! Transpose
      CD(i+1:Nseis,i)=0.d0                ! Clear lower triangle
    enddo
    deallocate(CDinv)
#endif
#ifdef REAL4
    CD4=real(CD,4)
#endif
    ! Write to file for debugging
    if(debugMode>1)then
      ! Save L1 parameters and Inverse of Cholesky decomposition of CD (binary)
      open(358,FILE='tempCDinv.tmp',ACCESS='stream',FORM='unformatted')
      write(358) int(NRseis,4)
      write(358) int(nT,4)
      write(358) int(Nseis,4)
      write(358) int(L1_all(1:NRseis),4)
      write(358) real(dble(L1_all(1:NRseis))*dt,4)
      do i=1,Nseis
        write(358) real(CD(i,1:Nseis),4)
      enddo
      close(358)
    endif
    
    deallocate(L1_all)
    
!---------------------------------------------------------------------
END SUBROUTINE



!---------------------------------------------------------------------
END !InitCore()



SUBROUTINE axcf(N,f1,f2,L1,L12,C)
!---------------------------------------------------------------------
!    Compute AXCF of functions f1 and f2 of samples N
!    L1 and L12 are time-shift windows in samples
!    -> Output into covariance matrix C
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: N,L1,L12
    real(8),intent(in):: f1(N),f2(N)
    real(8),intent(out):: C(N,N)
    complex(8),allocatable:: smooth1(:),smooth12(:),s1(:),s2(:),s1smooth1(:),s2smooth1(:),s2smooth12(:),XC(:,:)
    integer:: i,j,FFT,Lzeros
    
! Check the length of smoothing window L1
    if(L1<2)then
      write(*,*) ' Warning! L1 is too short'
    endif
    
! Allocation and init of variables
    FFT = 2**(int(log(dble(2*N))/log(2.d0)+0.99999999d0)+1)
    allocate(smooth1(FFT),s1(FFT),s2(FFT),s1smooth1(FFT),s2smooth1(FFT),s2smooth12(FFT),XC(FFT,FFT))
    C=0.d0
    s1=0.
    s2=0.
    s1(1:N) = f1(1:N)
    s2(1:N) = f2(1:N)
    
! Put zeros (values) before and after the signal
    Lzeros = max(L1,L12);
    s1(N+1:N+Lzeros) = f1(N)
    s2(N+1:N+Lzeros) = f2(N)
    s1(FFT-Lzeros+1:FFT) = f1(1)
    s2(FFT-Lzeros+1:FFT) = f2(1)
    
! f2 smoothing by L12 window
    if(L12>1)then
      allocate(smooth12(FFT))
      smooth12=0.
      smooth12(1:int(dble(L12)/2.+.51))=1.d0/dble(L12*FFT)
      smooth12(FFT:FFT-L12/2+1:-1)=1.d0/dble(L12*FFT)
      call four1(smooth12,FFT,1)
      call four1(s2,FFT,1)
      s2smooth12=s2*smooth12
      call four1(s2smooth12,FFT,-1)
      deallocate(smooth12)
    else
      s2smooth12=s2
    endif
    
! smoothing window L1
    smooth1=0.
    smooth1(1:int(dble(L1)/2.+.51))=1.d0/dble(L1*FFT)
    smooth1(FFT:FFT-L1/2+1:-1)=1.d0/dble(L1*FFT)
    call four1(smooth1,FFT,1)
    
! compute XCF
    do i=1,FFT
      XC(:,i)=s1(:)*cshift(s2smooth12(:),SHIFT=(FFT/2+1)-i)
      call four1(XC(:,i),FFT,1)
      XC(:,i)=XC(:,i)*smooth1(:)
      call four1(XC(:,i),FFT,-1)
    enddo

    call four1(s1,FFT,1)
    s1smooth1=s1*smooth1
    call four1(s1smooth1,FFT,-1)
    call four1(s2smooth12,FFT,1)
    s2smooth1=s2smooth12*smooth1
    call four1(s2smooth1,FFT,-1)
    
    do i=1,FFT
      XC(:,i) = XC(:,i) - s1smooth1(:)*cshift(s2smooth1(:),SHIFT=(FFT/2+1)-i)
    enddo
    
! Fill the covariance matrix by AXCF
    do i=1,N
     do j=1,N
       C(i,j)=dble(XC(i,FFT/2+1-(j-i)))
     enddo
    enddo

    deallocate(smooth1,s1,s2,s2smooth12,s1smooth1,s2smooth1,XC)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE saxcf(N,f1,f2,L1,L12,T,C)
!---------------------------------------------------------------------
!    Compute SAXCF of functions f1 and f2 of samples N
!    L1 and L12 are time-shift windows in samples
!    T is characteristic length of the signal in samples
!    -> Output into covariance matrix C
!---------------------------------------------------------------------
    implicit none
    real(8),parameter:: PI=3.1415926535d0
    integer,intent(in):: N,L1,L12,T
    real(8),intent(in):: f1(N),f2(N)
    real(8),intent(out):: C(N,N)
    complex(8),allocatable:: s1(:),s2(:),Rfg(:),fSACF(:),smooth1(:),smooth12(:)
    real(8),allocatable:: SACF(:),tw(:)
    real(8):: taper,r
    integer:: i,FFT

! Check the length of smoothing window L1
    if(L1<2)then
      write(*,*) ' Warning! L1 is too short'
    endif
    
! Allocation and init of variables
    FFT = 2**(int(log(dble(2*N))/log(2.d0)+0.99999999d0)+1)
    allocate(s1(FFT),s2(FFT),Rfg(FFT),fSACF(FFT),smooth1(FFT))
    allocate(SACF(2*N-1),tw(2*N-1))
    C=0.d0
    s1=0.
    s2=0.
    s1(1:N) = f1(1:N)
    s2(1:N) = f2(1:N)
    
! FFT of the f and g functions
    call four1(s1,FFT,1)
    call four1(s2,FFT,1)
    
! Prepare smoothing functions in freq. domain
    smooth1=0.
    smooth1(1:int(dble(L1)/2.+.51))=1.d0/dble(L1)
    smooth1(FFT:FFT-L1/2+1:-1)=1.d0/dble(L1)
    call four1(smooth1,FFT,1)
    
!  Cross-correlation of signals in freq. domain
    Rfg = conjg(s1)*s2
    
!  Convolution of triangle function (width 2*L1) and cross-correlation in freq. domain
     fSACF = Rfg - smooth1*conjg(smooth1)*Rfg

!  Smoothing by joint shift (width L12) in freq. domain
    if(L12>1)then
      allocate(smooth12(FFT))
      smooth12=0.
      smooth12(1:int(dble(L12)/2.+.51))=1.d0/dble(L12)
      smooth12(FFT:FFT-L12/2+1:-1)=1.d0/dble(L12)
      call four1(smooth12,FFT,1)
      fSACF = smooth12*fSACF
      deallocate(smooth12)
    endif
    
!  Norm by the effective signal length
    fSACF = fSACF/dble(T)
    
!  Back to time domain
    call four1(fSACF,FFT,-1)
    fSACF = fSACF/FFT
    SACF(1:N-1) = dble(fSACF(FFT-N+2:FFT))
    SACF(N:2*N-1) = dble(fSACF(1:N))
    
!  Tapered cosine window (SACF)
    taper = 0.666d0
    tw=1.d0
    do i=1,2*N-1
      r = dble(i-1)/dble(2*N-2)
      if(r.le.(taper/2.d0))then
        tw(i) = 0.5d0*( 1 + cos((2.d0*PI/taper)*(r-taper/2.d0)) )
      elseif(r.ge.(1.d0-taper/2.d0))then
        tw(i) = 0.5d0*( 1 + cos((2.d0*PI/taper)*(r-1+taper/2.d0)) )
      endif
    enddo
    tw = tw**3 ! cube for bigger effect
    SACF = SACF * tw !  Taper cosine window
    
! Fill the covariance matrix by SACF
    do i=0,N-1
      C(i+1,1:N) = SACF(N-i:2*N-i-1)
    enddo
    
    deallocate(s1,s2,Rfg,fSACF,SACF,tw,smooth1)

!---------------------------------------------------------------------
END SUBROUTINE


