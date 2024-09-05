!---------------------------------------------------------------------
!
!  Execute PSI inversion
!
!---------------------------------------------------------------------
!
!  Author: Miroslav Hallo
!  Charles University, Faculty of Mathematics and Physics
!  E-mail: hallo@karel.troja.mff.cuni.cz
!  Revision 10/2017: The initial version
!  Revision 2018: Fixed Yoffe error, improved random walker, hypocentrum
!  Revision 2019: Added REAL4, corrected reciprocal prior in birth/death
!
!  Copyright (C) 2017-2019  Miroslav Hallo
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

MODULE ExeBox
!---------------------------------------------------------------------
!  Module containing global model parameters
!---------------------------------------------------------------------
    ! Global model parameters (last dim for MC-chains)
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:):: m_veloc,m_risetime,m_slip
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:):: m_peaktime,m_rakes
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: m_hypo
    INTEGER,ALLOCATABLE,DIMENSION(:,:):: m_Nspline
    INTEGER,ALLOCATABLE,DIMENSION(:):: splineN0
    ! Global slip model parameters
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: slip,risetime,ruptime
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: peaktime,rakes
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:):: timefunc
#ifdef REAL4        
    REAL(4),ALLOCATABLE:: Dvec(:)
#else
    REAL(8),ALLOCATABLE:: Dvec(:)
#endif
    REAL(8):: moment
    ! PT global control
    REAL(8),ALLOCATABLE:: logPPDstore(:),misfitstore(:)
    INTEGER:: pert,iseed,nproc,rank,ifile,afile
    INTEGER:: abins,aaccept(2),aerrors(3),athr(8),aNspli(3)
    REAL(8):: aVR,sqUvec
    REAL,EXTERNAL:: ran3,gasdev !Random number generators
!---------------------------------------------------------------------
END MODULE



SUBROUTINE InvExe()
!---------------------------------------------------------------------
!  Creates initial model parameters and execute inversion,
!  the model space is prospected by a Monte-Carlo method
!---------------------------------------------------------------------
    use InputBox
    use ExeBox
    implicit none
    integer:: i
    integer:: ialg,nsteps,nchains,iburn,iseed0,nbins
    real(8):: swaprate,tlow,thigh
    character(len=80):: dir
    real(8),allocatable:: Chaintemp(:) ! Temperatures of each chain
    
    !----------------------------------
    ! Set up PT control parameters
    ialg = 0            ! Algorithm type: 
                        ! ialg=0, Parallel Tempering with T swap between all levels 
                        ! ialg=2, Parallel Tempering with T swap only allowed between 
                        ! neighbouring temperature levels
    nchains = nchainsG  ! Define number of chains
    swaprate = 1.0d0    ! Rate at which exchange swaps are proposed relative to within
                        ! chain steps. 1.0 = one exchage swap proposed for every McMC step.
                        ! Set this value to zero to turn off Parallel Tempering altogether.
    nsteps   = nstepsG  ! Number of chain steps per temperature
    iburn    = iburnG   ! Number of burn in samples
    iseed0 = 61254557   ! Random number seed
    tlow  = 1.d0        ! Lowest temperature of chains
    thigh = 50.d0       ! Highest temperature of chains
    nbins = 10          ! Number of temperature bins for diagnostics
    dir = './inv/'      ! Home directory for I/O files (some systems want full path names under MPI)
    
    ! Set global variables
    abins = 10         ! Number of steps (as one bin) for accept rate diagnostics
    
    !----------------------------------
    ! Print some drive variables
    if(debugMode>0)then
      write(*,*) " Number of cell points (along strike, along dip):"
      do i=1,NSeg
        write(*,'(2I5)') NL(i),NW(i)
      enddo
      write(*,*) " Number of control points (along strike, along dip):"
      do i=1,NSeg
        write(*,'(2I5)') NPL(i),NPW(i)
      enddo
      write(*,*) " Number of tigr points (along strike, along dip):"
      do i=1,NSeg
        write(*,'(2I6)') tigrNL(i),tigrNW(i)
      enddo
      write(*,*) " Spacing cell points [km] (along strike, along dip):"
      do i=1,NSeg
        write(*,'(2F9.3)') dL(i)/1000.d0,dW(i)/1000.d0
      enddo
      write(*,*) " Spacing control points [km] (along strike, along dip):"
      do i=1,NSeg
        write(*,'(2F9.3)') dPL(i)/1000.d0,dPW(i)/1000.d0
      enddo
      write(*,*) " Tigr spacing [km]:"
      write(*,'(F9.3)') tigrDD/1000.
      write(*,*) " Segment time delay [sec]:"
      do i=1,NSeg
        write(*,'(F9.3)') segDelay(i)
      enddo
      write(*,*) " Max. number of efficient spline control points:"
      write(*,'(I6)') MaxNSpline
    endif
    
    !----------------------------------
    ! Allocate
    allocate(Chaintemp(nchains))             ! Array for temperatures of each chain
    
    !----------------------------------
    ! Initialize PT library and determine process ID (if parallel)
    ! If we are in serial mode, then rank=0,nproc=1
    call PT (0,0,0,0,0,0.D0,0.D0,0.D0,0,0.0,0,dir,nproc,rank)
    
    !----------------------------------
    ! Initial model parameters
    write(*,*)' Initial model parameters'
    call InitMod(nchains,dir)
    
    !----------------------------------
    ! Print some PT control parameters
    if(debugMode>0)then
      write(*,*) " Number of cores and core rank:"
      write(*,'(2I5)') nproc,rank
      write(*,*) " Parallel tempering output directory:"
      write(*,'(A80)') dir
      write(*,*) " Limit temperatures of chains:"
      write(*,'(2F7.1)') tlow,thigh
      write(*,*) " Number of chains on different temperatures:"
      write(*,'(I7)') nchains
      write(*,*) " Rate at which exchange swaps within steps:"
      write(*,'(F8.2)') swaprate
      write(*,*) " Number of steps and burn steps:"
      write(*,'(2I8)') nsteps,iburn
    endif
    
    !----------------------------------
    ! Set temperatures for each chain (array Chaintemp)
    ! Random T-ladder with a log-uniform distribution between tlow and thigh
    write(*,*)' Set up temperatures'
    call Setuptempladder(nchains,tlow,thigh,dir,Chaintemp)
    
    !----------------------------------
    ! Calculate transition acceptance rates between chains
    ! at different temperatures and within chains at the same temperature
    write(*,*)' Parallel Tempering'
    call PT (1,ialg,nchains,nsteps,iburn,Chaintemp,thigh,tlow,&
              nbins,swaprate,iseed0,dir,nproc,rank)
    
    !----------------------------------
    ! Finish MPI and clean up
    write(*,*)' Finish and clean up'
    call PT (99,0,0,0,0,0.D0,0.D0,0.D0,0,0.0,0,dir,nproc,rank)
    
    close(ifile)
    close(afile)
    deallocate(Chaintemp)

!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE InitMod(nchains,dir)
!---------------------------------------------------------------------
!  Initial model parameters
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    implicit none
    integer,intent(in):: nchains
    character(len=80),intent(in):: dir
    character(len=4):: cproc
    character(len=96):: filename
    real(8):: a,sig
    integer im,kk,error,outThr,i,j
    real(8),allocatable:: tmp_moment(:)
    
    !----------------------------------
    ! Number of zero spline control points on edges
    allocate(splineN0(NSeg))
    do kk=1,NSeg
      if(ZSCtop==0)then
        splineN0(kk)=2*NL(kk)+2*NW(kk)
      else
        splineN0(kk)=NL(kk)+2*NW(kk)
      endif
    enddo
    
    !----------------------------------
    ! Initialize random number generator
    iseed=-(181615141+1000*rank*rank)
    
    !----------------------------------
    ! Alocate slip model arrays
    allocate(slip(maxval(NL),maxval(NW),NSeg))
    allocate(risetime(maxval(NL),maxval(NW),NSeg))
    allocate(ruptime(maxval(NL),maxval(NW),NSeg))
    allocate(peaktime(maxval(NL),maxval(NW),NSeg))
    allocate(rakes(maxval(NL),maxval(NW),NSeg))
    allocate(timefunc(Ssvd,maxval(NL),maxval(NW),NSeg))
    allocate(Dvec(Nseis))
    
    !----------------------------------
    ! Zeros into allocated slip model
    slip=0.d0
    risetime=0.d0
    ruptime=0.d0
    peaktime=0.d0
    rakes=0.d0
    timefunc=0.d0
    Dvec=0.d0
    
    !----------------------------------
    ! Allocate model arrays
    ! The last chain is artificial (nchains+1) used for perturbed model in cycles
    pert = nchains+1
    allocate(m_veloc(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_risetime(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_peaktime(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_rakes(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_hypo(3,NSeg,pert))
    allocate(m_slip(maxval(splineN0)+MaxNSpline,3,NSeg,pert))
    allocate(m_Nspline(NSeg,pert))
    ! Allocate model temp arrays
    allocate(logPPDstore(nchains))
    allocate(misfitstore(nchains))
    allocate(tmp_moment(nchains))
    
    !----------------------------------
    ! Spline zero control points on edges
    m_slip=0.d0 ! Required for zero control points on edges
    do im=1,nchains+1 ! Set for all including perturbed model
      do kk=1,NSeg
        if(ZSCtop==0)then
          m_slip(1:NL(kk),1,kk,im)=cellL(1:NL(kk),kk)
          m_slip(NL(kk)+1:2*NL(kk),1,kk,im)=cellL(1:NL(kk),kk)
          m_slip(NL(kk)+1:2*NL(kk),2,kk,im)=widt(kk)
          m_slip(2*NL(kk)+1:2*NL(kk)+NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
          m_slip(2*NL(kk)+NW(kk)+1:2*NL(kk)+2*NW(kk),1,kk,im)=leng(kk)
          m_slip(2*NL(kk)+NW(kk)+1:2*NL(kk)+2*NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
        else
          m_slip(1:NL(kk),1,kk,im)=cellL(1:NL(kk),kk)
          m_slip(NL(kk)+1:NL(kk)+NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
          m_slip(NL(kk)+NW(kk)+1:NL(kk)+2*NW(kk),1,kk,im)=leng(kk)
          m_slip(NL(kk)+NW(kk)+1:NL(kk)+2*NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
        endif
      enddo
    enddo
    
    !----------------------------------
    ! Open output file for diagnostics
    write(cproc,'(I0)') rank
    filename=trim(dir)//'log/init_'//trim(cproc)//'.log'
    open(unit=48,file=filename,action='write',status='replace')
    write(48,'(A)') 'Initial models for all chains'
    write(48,'(A)') 'Initial parameters veloc, peaktime, risetime, rakes are constant over segment'
    write(48,'(A,I8)') 'Total number of chains:',nchains
    write(48,'(A)') '----------------------------------------'
    
    !----------------------------------
    ! Initial random model parameters
    do im=1,nchains
      do while(.true.) ! until no error
        outThr = 0
        do kk=1,NSeg
        
          ! Init random velocity models (gauss)
          sig = (thr_velo(2)-thr_velo(1))/4.d0
          a = (thr_velo(1) + (thr_velo(2)-thr_velo(1))/2.d0) + dble(gasdev(iseed)) * sig
          if (a<thr_velo(1))then
            a=thr_velo(1)
            outThr = 1
          elseif(a>thr_velo(2))then
            a=thr_velo(2)
            outThr = 1
          endif
          m_veloc(1:NPL(kk),1:NPW(kk),kk,im) = a
        
          ! Init random risetime (uniform)
          a = thr_rise(1) + ran3(iseed)*(thr_rise(2)-thr_rise(1))
          m_risetime(1:NPL(kk),1:NPW(kk),kk,im) = a
        
          ! Init random peaktime (uniform)
          a = thr_peak(1) + ran3(iseed)*(thr_peak(2)-thr_peak(1))
          m_peaktime(1:NPL(kk),1:NPW(kk),kk,im) = a
          
          ! Check for Yoffe function error (risetime < 1.575*peaktime)
          if ( m_risetime(1,1,kk,im) < (1.575d0*m_peaktime(1,1,kk,im) + 3.d0*abs(thr_peak(3))) )then
            outThr = 3
          endif
          
          ! Init random rakes (gauss)
          if(nrake==1)then
            a = rake(kk)
          elseif(nrake==2)then
            sig = thr_rake(1)/2.d0
            a = rake(kk) + dble(gasdev(iseed)) * sig
            if(a<rake(kk)-thr_rake(1))then
              a=rake(kk)-thr_rake(1)
              outThr = 4
            elseif(a>rake(kk)+thr_rake(1))then
              a=rake(kk)+thr_rake(1)
              outThr = 4
            endif
            if(a<rake(kk).and.a<-180.d0)then
              a=a+360.d0
            elseif(a>rake(kk).and.a>180.d0)then
              a=a-360.d0
            endif
          endif
          m_rakes(1:NPL(kk),1:NPW(kk),kk,im) = a
        
          ! Init random hypocentrum-leng (gauss)
          sig = thr_hypo(1)/2.d0
          a = epicL(kk) + dble(gasdev(iseed)) * sig
          if (a<epicL(kk)-thr_hypo(1))then
            a=epicL(kk)-thr_hypo(1)
            outThr = 5
          elseif(a>epicL(kk)+thr_hypo(1))then
            a=epicL(kk)+thr_hypo(1)
            outThr = 5
          endif
          if (a<0.d0)then
            a=0.d0
            outThr = 5
          elseif(a>leng(kk))then
            a=leng(kk)
            outThr = 5
          endif
          m_hypo(1,kk,im) = a
          
          ! Init random hypocentrum-widt (gauss)
          sig = thr_hypo(2)/2.d0
          a = epicW(kk) + dble(gasdev(iseed)) * sig
          if (a<epicW(kk)-thr_hypo(2))then
            a=epicW(kk)-thr_hypo(2)
            outThr = 5
          elseif(a>epicW(kk)+thr_hypo(2))then
            a=epicW(kk)+thr_hypo(2)
            outThr = 5
          endif
          if (a<0.d0)then
            a=0.d0
            outThr = 5
          elseif(a>widt(kk))then
            a=widt(kk)
            outThr = 5
          endif
          m_hypo(2,kk,im) = a
          
          ! Init random hypocentrum-delay (uniform)
          if(kk.eq.1)then
            a = segDelay(1) + (ran3(iseed)-0.5)*thr_delay(1)*2.d0
          else
            a = segDelay(kk) + (ran3(iseed)-0.5)*thr_delay(2)*2.d0
          endif
          if (a<0.d0)then
            a=0.d0
            outThr = 6
          endif
          m_hypo(3,kk,im) = a
          
          ! Init splines number (static)
          m_Nspline(kk,im) = 1
          
          ! Init random spline1-leng (gauss)
          sig = (leng(kk)-dL(kk)) / 4.d0
          a = (leng(kk)/2.d0) + dble(gasdev(iseed)) * sig
          if (a<(dL(kk)/2.d0))then
            a=dL(kk)/2.d0
            outThr = 7
          elseif(a>(leng(kk)-(dL(kk)/2.d0)))then
            a=leng(kk)-(dL(kk)/2.d0)
            outThr = 7
          endif
          m_slip(splineN0(kk)+1,1,kk,im) = a
          
          ! Init random spline1-widt (gauss)
          sig = (widt(kk)-dW(kk)) / 4.d0
          a = (widt(kk)/2.d0) + dble(gasdev(iseed)) * sig
          if (a<(dW(kk)/2.d0))then
            a=(dW(kk)/2.d0)
            outThr = 7
          elseif(a>(widt(kk)-(dW(kk)/2.d0)))then
            a=widt(kk)-(dW(kk)/2.d0)
            outThr = 7
          endif
          m_slip(splineN0(kk)+1,2,kk,im) = a
          
          ! Init spline1-slip (static)
          m_slip(splineN0(kk)+1,3,kk,im) = 1.d0
          
        enddo
        
        ! Init random seismic moment (gauss)
        sig = thr_M0/2.d0
        a = Mfix + dble(gasdev(iseed)) * sig
        if (a<(Mfix-thr_M0))then
          a=Mfix-thr_M0
          outThr = 8
        elseif(a>(Mfix+thr_M0))then
          a=Mfix+thr_M0
          outThr = 8
        endif
        tmp_moment(im) = a
        
        ! --------------------------
        ! Load initial model from file (for debugging and so..)
        ! (discard all generated initial parameters)
        if(.false.)then
          filename='inv/best0_model.txt'
          open(237,FILE=filename)
          do kk=1,NSeg
            read(237,*)
            read(237,*)
            do j=1,NPW(kk)
              read(237,*) (m_veloc(i,j,kk,im) ,i=1,NPL(kk))
            enddo
            read(237,*)
            do j=1,NPW(kk)
              read(237,*) (m_risetime(i,j,kk,im) ,i=1,NPL(kk))
            enddo
            read(237,*)
            do j=1,NPW(kk)
              read(237,*) (m_peaktime(i,j,kk,im) ,i=1,NPL(kk))
            enddo
            read(237,*)
            do j=1,NPW(kk)
              read(237,*) (m_rakes(i,j,kk,im) ,i=1,NPL(kk))
            enddo
            read(237,*)
            read(237,*) m_hypo(1,kk,im),m_hypo(2,kk,im),m_hypo(3,kk,im)
            read(237,*)
            read(237,*) j, m_Nspline(kk,im)
            read(237,*)
            read(237,*) (m_slip(i,1,kk,im) ,i=1,j)
            read(237,*) (m_slip(i,2,kk,im) ,i=1,j)
            read(237,*) (m_slip(i,3,kk,im) ,i=1,j)
          enddo
          close(237)
          tmp_moment(im) = Mfix
          outThr = 0
        endif
        
        ! Cycle if any model parameter on the edge
        if(outThr.ne.0)then
          write(48,'(A9,I1,A11,I8)') 'GenThr[',outThr,'] in chain:',im
          flush(48)
          cycle
        endif
        
        ! Calculate -logPPD of initial model (with unit slip)
        call Forward(im,tmp_moment(im),misfitstore(im),logPPDstore(im),error)
        
        ! Cycle if any error
        if(error.ne.0)then
          write(48,'(A9,I1,A11,I8)') 'GenErr[',error,'] in chain:',im
          flush(48)
          cycle
        endif
        
        ! If no error
        if(error.eq.0.and.outThr.eq.0)then
          exit
        endif
        
      enddo ! until no error
    enddo
    
    ! Prepare vars for diagnostics
    aaccept = 0
    aerrors = 0
    athr = 0
    aNspli = 0
    aVR = 0.d0
    sqUvec = dot_product(Uvec,Uvec)
    
    ! Write resultant initial models for diagnostics
    do kk=1,NSeg
      write(48,'(A)') '----------------------------------------'
      write(48,'(A13,I3)') 'Segment No.',kk
      write(48,'(A18,100000F10.1)') 'Velocity [m/s]:',m_veloc(1,1,kk,1:nchains)
      write(48,'(A18,100000F10.2)') 'Risetime [s]:',m_risetime(1,1,kk,1:nchains)
      write(48,'(A18,100000F10.2)') 'Peaktime [s]:', m_peaktime(1,1,kk,1:nchains)
      write(48,'(A18,100000F10.1)') 'Rakes [deg]:',m_rakes(1,1,kk,1:nchains)
      write(48,'(A18,100000F10.1)') 'Hypo-leng [km]:',m_hypo(1,kk,1:nchains)/1000.d0
      write(48,'(A18,100000F10.1)') 'Hypo-widt [km]:',m_hypo(2,kk,1:nchains)/1000.d0
      write(48,'(A18,100000F10.2)') 'Hypo-delay [s]:',m_hypo(3,kk,1:nchains)
      write(48,'(A18,100000I10)') 'Splines No.:',m_Nspline(kk,1:nchains)
      write(48,'(A18,100000F10.1)') 'Spline1-leng [km]:',m_slip(splineN0(kk)+1,1,kk,1:nchains)/1000.d0
      write(48,'(A18,100000F10.1)') 'Spline1-widt [km]:',m_slip(splineN0(kk)+1,2,kk,1:nchains)/1000.d0
      write(48,'(A18,100000F10.3)') 'Spline1-val [m]:',m_slip(splineN0(kk)+1,3,kk,1:nchains)
    enddo
    write(48,'(A)') '----------------------------------------'
    write(48,'(A7)') 'Total'
    write(48,'(A18,100000E10.3)') 'Seis.moment [N-m]:',tmp_moment(1:nchains)
    write(48,'(A18,100000F10.3)') 'Mw:',(log10(tmp_moment(1:nchains))-9.1d0)/1.5d0
    write(48,'(A18,100000F10.3)') 'VR [%]:',(1-(misfitstore(1:nchains)/sqUvec))*100.d0
    close(48)
    
    !----------------------------------
    ! Open output file for models
    ifile = 55
    filename=trim(dir)//'xmodels_'//trim(cproc)//'.dat'
    open(unit=ifile,form='unformatted',file=filename,action='write',status='replace')
    flush(ifile)
    
    !----------------------------------
    ! Open output file for diagnostics
    afile = 69
    filename=trim(dir)//'log/stats_'//trim(cproc)//'.log'
    open(unit=afile,file=filename,action='write',status='replace')
    write(afile,'(A)') 'Run-time statistics from parallel tempering'
    write(afile,'(A)') 'Legend:'
    write(afile,'(A)') '  AcceptT1 - Number of accepted new models on T=1'
    write(afile,'(A)') '  AcceptAll - Number of accepted new models on all Ts'
    write(afile,'(A)') '  GenErr[1] - Rupture-time generation Err (time_2D)'
    write(afile,'(A)') '  GenErr[2] - Spline inverse Err (dgesv)'
    write(afile,'(A)') '  GenErr[3] - Yoffe function Err (risetime < 1.575*peaktime)'
    write(afile,'(A)') '  VR[%] - The best Variance reduction from all chains'
    write(afile,'(A)') '  GenThr[1] - Threshold reached in velocity'
    write(afile,'(A)') '  GenThr[2] - Threshold reached in risetime'
    write(afile,'(A)') '  GenThr[3] - Threshold reached in peaktime'
    write(afile,'(A)') '  GenThr[4] - Threshold reached in rakes'
    write(afile,'(A)') '  GenThr[5] - Threshold reached in hypo-pos'
    write(afile,'(A)') '  GenThr[6] - Threshold reached in hypo-delay'
    write(afile,'(A)') '  GenThr[7] - Threshold reached in spline-pos/val'
    write(afile,'(A)') '  GenThr[8] - Info: The seismic moment exceeded the assumed value'
    write(afile,'(A)') '------------------------------------------------------------'
    write(afile,'(18A10)') 'CHsteps','AcceptT1','AcceptAll','GenErr[1]','GenErr[2]','GenErr[3]', &
                           'VR[%]','SplinAll','SplinKill','SplinNew', &
             'GenThr[1]','GenThr[2]','GenThr[3]','GenThr[4]','GenThr[5]','GenThr[6]','GenThr[7]','GenThr[8]'
    write(afile,'(A10,5I10,12A10)') 'Max:',abins,abins*nchains,abins*nchains,abins*nchains,abins*nchains, &
                  '100%','per_segm.','propos %','propos %','%','%','%','%','%','%','%','%'
                  
    write(afile,'(A)') '------------------------------------------------------------'
    flush(afile)
    
    deallocate(tmp_moment)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE AdvanceChain(ichain,is,Temper,logPPD)
!---------------------------------------------------------------------
!  User routine to advance a Metropolis-Hastings random 
!  walker using Markov chain Monte Carlo.
!  In:   ichain     Integer   : Random walk (chain) identifer
!        is         Integer   : Main loop ID over chain steps (iburn+nsteps)
!        Temper     Doubler   : Temperature of chain supplied by the calling routine
!  Out:  logPPD     Real*8    : Negative logrithm of the posterior
!                               probability density function of updated chain.
!---------------------------------------------------------------------
    use InputBox
    use ExeBox
    use LibBox
    implicit none
    integer,intent(in):: ichain,is
    real(8),intent(in):: Temper
    real(8),intent(out):: logPPD
    
    integer:: i,j,kk,error
    real(8):: misfit,a,sig,thr,logPPD1,logPPD2
    real(8):: logQ12,logQ21,q12,q21,tmp
    logical:: yn
    character(len=19):: filename
    
    !----------------------------------
    !  Get negative log PPD of current model
    logPPD1 = logPPDstore(ichain)
    logPPD = logPPD1
    
    !----------------------------------
    !  Calculate perturbed model as +-1 perturbation from existing model
    do kk=1,NSeg
      do j=1,NPW(kk)
        do i=1,NPL(kk)
        
          ! Random velocity step (gauss)
          a = m_veloc(i,j,kk,ichain) + dble(gasdev(iseed)) * thr_velo(3)
          if (a<thr_velo(1))then
            a=(2.d0*thr_velo(1))-a
            if(a>thr_velo(2))then ! For sure
              a=m_veloc(i,j,kk,ichain)
            endif
            athr(1)=athr(1)+1
          elseif(a>thr_velo(2))then
            a=(2.d0*thr_velo(2))-a
            if(a<thr_velo(1))then ! For sure
              a=m_veloc(i,j,kk,ichain)
            endif
            athr(1)=athr(1)+1
          endif
          m_veloc(i,j,kk,pert) = a
          
          ! Random risetime step (gauss)
          a = m_risetime(i,j,kk,ichain) + dble(gasdev(iseed)) * thr_rise(3)
          if (a<thr_rise(1))then
            a=(2.d0*thr_rise(1))-a
            if(a>thr_rise(2))then ! For sure
              a=m_risetime(i,j,kk,ichain)
            endif
            athr(2)=athr(2)+1
          elseif(a>thr_rise(2))then
            a=(2.d0*thr_rise(2))-a
            if(a<thr_rise(1))then ! For sure
              a=m_risetime(i,j,kk,ichain)
            endif
            athr(2)=athr(2)+1
          endif
          m_risetime(i,j,kk,pert) = a
          
          ! Random peaktime step (gauss)
          a = m_peaktime(i,j,kk,ichain) + dble(gasdev(iseed)) * thr_peak(3)
          if (a<thr_peak(1))then
            a=(2.d0*thr_peak(1))-a
            if(a>thr_peak(2))then ! For sure
              a=m_peaktime(i,j,kk,ichain)
            endif
            athr(3)=athr(3)+1
          elseif(a>thr_peak(2))then
            a=(2.d0*thr_peak(2))-a
            if(a<thr_peak(1))then ! For sure
              a=m_peaktime(i,j,kk,ichain)
            endif
            athr(3)=athr(3)+1
          endif
          m_peaktime(i,j,kk,pert) = a
          
          ! Random rake step (gauss)
          if(nrake==1)then
            a = m_rakes(i,j,kk,ichain)
          elseif(nrake==2)then
            a = m_rakes(i,j,kk,ichain) + dble(gasdev(iseed)) * thr_rake(2)
            sig = abs(rake(kk)-a)
            if(sig.gt.180.d0)then
              sig = abs(360.d0 - sig)
            endif
            if(sig.gt.thr_rake(1))then
              a=m_rakes(i,j,kk,ichain)
              athr(4)=athr(4)+1
            endif
            if(a<m_rakes(i,j,kk,ichain).and.a<-180.d0)then
              a=a+360.d0
            elseif(a>m_rakes(i,j,kk,ichain).and.a>180.d0)then
              a=a-360.d0
            endif
          endif
          m_rakes(i,j,kk,pert) = a
          
        enddo
      enddo
      
      ! Random hypocentrum-leng step (gauss)
      a = m_hypo(1,kk,ichain) + dble(gasdev(iseed)) * thr_hypo(3)
      if (a<epicL(kk)-thr_hypo(1))then
        a=(2.d0*(epicL(kk)-thr_hypo(1)))-a
        if(a>epicL(kk)+thr_hypo(1))then ! For sure
          a=m_hypo(1,kk,ichain)
        endif
        athr(5)=athr(5)+1
      elseif(a>epicL(kk)+thr_hypo(1))then
        a=(2.d0*(epicL(kk)+thr_hypo(1)))-a
        if(a<epicL(kk)-thr_hypo(1))then ! For sure
          a=m_hypo(1,kk,ichain)
        endif
        athr(5)=athr(5)+1
      endif
      if (a<0.d0)then
        a=abs(a)
        if(a>epicL(kk)+thr_hypo(1))then ! For sure
          a=m_hypo(1,kk,ichain)
        elseif(a>leng(kk))then
          a=m_hypo(1,kk,ichain)
        endif
        athr(5)=athr(5)+1
      elseif(a>leng(kk))then
        a=(2.d0*leng(kk))-a
        if(a<epicL(kk)-thr_hypo(1))then ! For sure
          a=m_hypo(1,kk,ichain)
        elseif(a<0.d0)then
          a=m_hypo(1,kk,ichain)
        endif
        athr(5)=athr(5)+1
      endif
      m_hypo(1,kk,pert) = a
      
      ! Random hypocentrum-widt step (gauss)
      a = m_hypo(2,kk,ichain) + dble(gasdev(iseed)) * thr_hypo(3)
      if (a<epicW(kk)-thr_hypo(2))then
        a=(2.d0*(epicW(kk)-thr_hypo(2)))-a
        if(a>epicW(kk)+thr_hypo(2))then ! For sure
          a=m_hypo(2,kk,ichain)
        endif
        athr(5)=athr(5)+1
      elseif(a>epicW(kk)+thr_hypo(2))then
        a=(2.d0*(epicW(kk)+thr_hypo(2)))-a
        if(a<epicW(kk)-thr_hypo(2))then ! For sure
          a=m_hypo(2,kk,ichain)
        endif
        athr(5)=athr(5)+1
      endif
      if (a<0.d0)then
        a=abs(a)
        if(a>epicW(kk)+thr_hypo(2))then ! For sure
          a=m_hypo(2,kk,ichain)
        elseif(a>widt(kk))then
          a=m_hypo(2,kk,ichain)
        endif
        athr(5)=athr(5)+1
      elseif(a>widt(kk))then
        a=(2.d0*widt(kk))-a
        if(a<epicW(kk)-thr_hypo(2))then ! For sure
          a=m_hypo(2,kk,ichain)
        elseif(a<0.d0)then
          a=m_hypo(2,kk,ichain)
        endif
        athr(5)=athr(5)+1
      endif
      m_hypo(2,kk,pert) = a
      
      ! Random hypocentrum-delay step (gauss)
      a = m_hypo(3,kk,ichain) + dble(gasdev(iseed)) * thr_delay(3)
      if(kk.eq.1)then
        thr=thr_delay(1)
      else
        thr=thr_delay(2)
      endif
      if(a<segDelay(kk)-thr)then
        a=(2.d0*(segDelay(kk)-thr))-a
        if(a>segDelay(kk)+thr)then ! For sure
          a=m_hypo(3,kk,ichain)
        endif
        athr(6)=athr(6)+1
      elseif(a>segDelay(kk)+thr)then
        a=(2.d0*(segDelay(kk)+thr))-a
        if(a<segDelay(kk)-thr)then ! For sure
          a=m_hypo(3,kk,ichain)
        endif
        athr(6)=athr(6)+1
      endif
      if (a<0.d0)then
        a=abs(a)
        if(a>segDelay(kk)+thr)then ! For sure
          a=m_hypo(3,kk,ichain)
        endif
        athr(6)=athr(6)+1
      endif
      m_hypo(3,kk,pert) = a
      
      ! Existing splines
      m_Nspline(kk,pert) = m_Nspline(kk,ichain)
      do j=1,m_Nspline(kk,pert)
        aNspli(1) = aNspli(1) + 1
        
        ! Random splines-leng step (gauss)
        a = m_slip(splineN0(kk)+j,1,kk,ichain) + dble(gasdev(iseed)) * thr_spline(1)
        if (a<(dL(kk)/2.d0))then
          a=dL(kk)-a
          if(a>(leng(kk)-(dL(kk)/2.d0)))then ! For sure
            a=m_slip(splineN0(kk)+j,1,kk,ichain)
          endif
          athr(7)=athr(7)+1
        elseif(a>(leng(kk)-(dL(kk)/2.d0)))then
          a=(2.d0*(leng(kk)-(dL(kk)/2.d0)))-a
          if(a<(dL(kk)/2.d0))then ! For sure
            a=m_slip(splineN0(kk)+j,1,kk,ichain)
          endif
          athr(7)=athr(7)+1
        endif
        m_slip(splineN0(kk)+j,1,kk,pert) = a
      
        ! Random splines-widt step (gauss)
        a = m_slip(splineN0(kk)+j,2,kk,ichain) + dble(gasdev(iseed)) * thr_spline(1)
        if (a<(dW(kk)/2.d0))then
          a=dW(kk)-a
          if(a>(widt(kk)-(dW(kk)/2.d0)))then ! For sure
            a=m_slip(splineN0(kk)+j,2,kk,ichain)
          endif
          athr(7)=athr(7)+1
        elseif(a>(widt(kk)-(dW(kk)/2.d0)))then
          a=(2.d0*(widt(kk)-(dW(kk)/2.d0)))-a
          if(a<(dW(kk)/2.d0))then ! For sure
            a=m_slip(splineN0(kk)+j,2,kk,ichain)
          endif
          athr(7)=athr(7)+1
        endif
        m_slip(splineN0(kk)+j,2,kk,pert) = a
        
        ! Random splines-val step (gauss)
        a = m_slip(splineN0(kk)+j,3,kk,ichain) + dble(gasdev(iseed)) * thr_spline(2)
        !if (a<0.d0)then
        !  a=0.d0
        !  athr(7)=athr(7)+1
        !endif
        m_slip(splineN0(kk)+j,3,kk,pert) = a
        
      enddo
    enddo
    
    !----------------------------------
    ! Change dimension (trans-D)
    
    q12 = 1.d0  ! q(m1|m2)
    q21 = 1.d0  ! q(m2|m1)
    
    a = ran3(iseed)
    if(a.le.ptrans)then
      
      ! Select segment
      kk = 1 + int(ran3(iseed)*NSeg)
      if(kk.gt.NSeg)then
        kk=NSeg
      endif
    
      a = ran3(iseed)
      if(a.lt.0.5d0)then
        !----------------------------------
        ! Dead move
        !----------------------------------
        if(m_Nspline(kk,pert).gt.1)then ! Is death possible?
          aNspli(2) = aNspli(2) + 1
          ! Select spline to kill, replace by last and kill
          j = 1 + int(ran3(iseed)*m_Nspline(kk,pert))
          if(j.lt.m_Nspline(kk,pert))then
            m_slip(splineN0(kk)+j,1:3,kk,pert) = m_slip(splineN0(kk)+m_Nspline(kk,pert),1:3,kk,pert)
          endif
          
          ! Proposal distribution
          q12 = dble(m_Nspline(kk,pert))            ! q(m1|m2)
          q21 = dble(m_Nspline(kk,pert)-1)          ! q(m2|m1)
          m_Nspline(kk,pert) = m_Nspline(kk,pert)-1
        endif
        
      elseif(m_Nspline(kk,pert).lt.MaxNSpline)then
        !----------------------------------
        ! Birth move
        !----------------------------------
        aNspli(3) = aNspli(3) + 1
        j = splineN0(kk)+m_Nspline(kk,pert)
        
        ! Init random spline1-leng (uniform)
        a = (dL(kk)/2.d0) + ran3(iseed)*(leng(kk)-dL(kk))
        m_slip(j+1,1,kk,pert) = a
        
        ! Init random spline1-widt (uniform)
        a = (dW(kk)/2.d0) + ran3(iseed)*(widt(kk)-dW(kk))
        m_slip(j+1,2,kk,pert) = a
        
        ! Random splines-val step (gauss) from init position
        sig = sum(bspline(j,1,1,m_slip(1:j,1,kk,pert),m_slip(1:j,2,kk,pert),m_slip(1:j,3,kk,pert),&
            m_slip(j+1,1,kk,pert),m_slip(j+1,2,kk,pert)))
        a = sig + dble(gasdev(iseed)) * thr_spline(2)
        !if (a<0.d0)then
        !  a=0.d0
        !endif
        m_slip(j+1,3,kk,pert) = a
        
        ! Proposal distribution
        q12 = dble(m_Nspline(kk,pert))            ! q(m1|m2)
        q21 = dble(m_Nspline(kk,pert)+1)          ! q(m2|m1)
        m_Nspline(kk,pert) = m_Nspline(kk,pert)+1
        
      endif
    endif
    
    !----------------------------------
    ! Calculate -logPPD of proposed model in ichain
    call Forward(pert,-12345.d0,misfit,logPPD2,error)
    
    !----------------------------------
    ! Check for total seismic moment
    if(error.eq.0)then
      if (moment<(Mfix-thr_M0))then
        !error = 4
        athr(8)=athr(8)+1
      elseif(moment>(Mfix+thr_M0))then
        !error = 4
        athr(8)=athr(8)+1
      endif
    endif
    
    !----------------------------------
    ! log of probability of moving state
    logQ21 = log(q21)
    logQ12 = log(q12)
    
    !----------------------------------
    ! Decide whether to accept perturbed model using M-H criterion.
    ! raised to the power of inverse temperature.
    if(error.eq.0)then
      call PT_McMC_accept(Temper,logPPD1,logQ12,logPPD2,logQ21,yn)
    else
      yn = .false.
    endif
    
    !----------------------------------
    ! We accept the new model
    if(yn)then
      m_veloc(:,:,:,ichain) = m_veloc(:,:,:,pert)
      m_risetime(:,:,:,ichain) = m_risetime(:,:,:,pert)
      m_peaktime(:,:,:,ichain) = m_peaktime(:,:,:,pert)
      m_rakes(:,:,:,ichain) = m_rakes(:,:,:,pert)
      m_hypo(:,:,ichain) = m_hypo(:,:,pert)
      m_Nspline(:,ichain) = m_Nspline(:,pert)
      m_slip(:,:,:,ichain) = m_slip(:,:,:,pert)
      logPPDstore(ichain) = logPPD2
      misfitstore(ichain) = misfit
      logPPD = logPPD2
    endif
    
    !----------------------------------
    ! Write out target model
    if(Temper.eq.1.d0 .and. is.gt.iburnG)then
      misfit = misfitstore(ichain)
      write(ifile) logPPD,misfit ! 2X real(8)
      write(ifile) m_veloc(:,:,:,ichain)
      write(ifile) m_risetime(:,:,:,ichain)
      write(ifile) m_peaktime(:,:,:,ichain)
      write(ifile) m_rakes(:,:,:,ichain)
      write(ifile) m_hypo(:,:,ichain)
      write(ifile) m_Nspline(:,ichain)
      write(ifile) m_slip(:,:,:,ichain)
      flush(ifile)
    endif
    
    !----------------------------------
    ! Write out diagnostics
    if(yn)then
      aaccept(2)=aaccept(2)+1
      if(Temper.eq.1)then
        aaccept(1)=aaccept(1)+1
      endif
      if(aVR.lt.(1-(misfit/sqUvec)))then
        aVR = 1-(misfit/sqUvec)
      endif
      
    elseif(error.ne.0 .and. error.lt.4)then
      aerrors(error)=aerrors(error)+1
    endif
    if(mod(is,abins).eq.0.and.ichain.eq.1)then
      tmp = dble(abins*nchainsG)
      write(afile,'(6I10,4F10.2,8F10.2)') is,aaccept,aerrors,aVR*100.d0,aNspli(1)/(tmp*NSeg),100.*aNspli(2:3)/tmp,&
                    100.*athr(1:4)/(tmp*sum(NPW*NPL)),100.*athr(5)/(2*tmp*NSeg),&
                    100.*athr(6)/(tmp*NSeg),100.*athr(7)/(3*aNspli(1)),100.*athr(8)/tmp
      flush(afile)
      aaccept=0
      aerrors=0
      athr=0
      aNspli=0
      aVR=0.d0
    endif
    
    !----------------------------------
    ! Write a time-slice from a chain
    if(.FALSE.)then
    !if(ichain.eq.1 .and. (rank.eq.1 .or. nproc.eq.1) )then
      !  Filename
      if(is<10)then
        write(filename,'(A14,I1,A4)')'inv/mov/000000',is,'.txt'
      elseif(is<100)then
        write(filename,'(A13,I2,A4)')'inv/mov/00000',is,'.txt'
      elseif(is<1000)then
        write(filename,'(A12,I3,A4)')'inv/mov/0000',is,'.txt'
      elseif(is<10000)then
        write(filename,'(A11,I4,A4)')'inv/mov/000',is,'.txt'
      elseif(is<100000)then
        write(filename,'(A10,I5,A4)')'inv/mov/00',is,'.txt'
      elseif(is<1000000)then
        write(filename,'(A9,I6,A4)')'inv/mov/0',is,'.txt'
      else
        write(filename,'(A8,I7,A4)')'inv/mov/',is,'.txt'
      endif
      !  Total Slip time-slice
      kk=1
      open(292,file=filename)
      do j=1,NW(kk)
        write(292,'(1000E16.7E3)')(sum(timefunc(:,i,j,kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),j,kk))*dt
      enddo
      write(292,'(1000E16.7E3)')(sum(timefunc(:,i,NW(kk),kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),NW(kk),kk))*dt
      close(292)
    endif
    
    
    RETURN
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE Forward(im,targetM0,misfit,logPDF,error)
!---------------------------------------------------------------------
!  Calculate forward problem evaluate misfit and -log(PDF)
!  In:  im - model identifer
!       targetM0 - target total seismic moment (scale values in m_slip)
!  Out: misfit - Misfit of data and synthetics
!       logPDF - Negative logrithm of the posterior PDF
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    implicit none
    integer,intent(in):: im
    real(8),intent(in):: targetM0
    real(8),intent(out):: logPDF
    real(8),intent(out):: misfit
    integer,intent(out):: error
    
    !----------------------------------
    ! Interpret model parameters
    !  (into global vars: slip,risetime,ruptime,peaktime,rakes,moment)
    call InterpMod(im,targetM0,error)
    
    !----------------------------------
    !  If error
    if(error.ne.0)then
      return
    endif
    
    !----------------------------------
    ! Evaluate misfit function (90% time of the forward problem)
#ifdef MKL
#ifdef REAL4
    call sgemv('N',Nseis,Msvd,1.,H,Nseis,tfin,1,0.,Dvec,1)
#else
    call dgemv('N',Nseis,Msvd,1.d0,H,Nseis,tfin,1,0.d0,Dvec,1)
#endif
#else
    Dvec=matmul(H,tfin)
#endif
    misfit=dot_product(Dvec-Uvec,Dvec-Uvec)
    !misfit=0.d0
    
    !----------------------------------
    ! Calculate the minus logPDF of input model
    ! logPDF evaluated up to an additive constant, because in PT
    ! the decision to accept depends only on differences logPDF
    !logPDF = -log( exp(-0.5d0*misfit) )
    logPDF = (0.5d0*misfit)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE InterpMod(im,targetM0,error)
!---------------------------------------------------------------------
!  Model parameters into slip,risetime,ruptime,peaktime,rakes
!  It generates also slip-rate time source functions (tfin)
!  It proceed model with identifer im (last index of arrays)
!  It may linearly scale values of spline points (m_slip) to reach the
!     targetM0 (total seismic moment), to disable this feature put targetM0 <= 0.d0
!  error 0 - NO_ERROR
!  error 1 - TIME_2D PROBLEM
!  error 2 - SPLINE DGESV PROBLEM
!  error 3 - RYOFFE FUNCTION PROBLEM
!  error 4 - TARGET SEISMIC MOMENT DIDNT REACHED
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    use LibBox
    implicit none
    integer,intent(in):: im
    real(8),intent(in):: targetM0
    integer,intent(out):: error
    integer*4,external:: time_2d
    integer*4:: iost
    real(4):: xg,zg,eps_init
    real(4),allocatable,dimension(:):: hs,t_rupt
    real(8),allocatable,dimension(:,:):: tigr_res
    real(8),allocatable,dimension(:,:,:):: m_tmp_rakes
    real(8),allocatable,dimension(:):: rake0
    integer:: i,j,k,kk,iter,iterN,resh(2),tmp,r,coun
    real(8):: rangd
    
    error=0
    !----------------------------------
    ! Compute rupture time in ti(me)gr(id)
    do kk=1,NSeg
      allocate(tigr_res(tigrNL(kk),tigrNW(kk)))
      allocate(hs(tigrNL(kk)*tigrNW(kk)),t_rupt(tigrNL(kk)*tigrNW(kk)))
      t_rupt=0.
      ! Interpolate control points to tigr
      do j=1,tigrNW(kk)
        do i=1,tigrNL(kk)
          hs(i+(j-1)*tigrNL(kk))=tigrDD/real(interp2D(NPL(kk),cpL(1:NPL(kk),kk),NPW(kk),&
          cpW(1:NPW(kk),kk),m_veloc(1:NPL(kk),1:NPW(kk),kk,im),tigrL(i,kk),tigrW(j,kk)))
        enddo
      enddo
      ! Set param for eikonal solver
      eps_init=0.
      xg=real(m_hypo(1,kk,im))/tigrDD+1.
      zg=real(m_hypo(2,kk,im))/tigrDD+1.
      ! Run C++ Eikonal solver
      iost=time_2d(hs,t_rupt,tigrNL(kk),tigrNW(kk),xg,zg,eps_init,0)
      ! Interpolate tigr to cell grid
      resh(1)=tigrNL(kk)
      resh(2)=tigrNW(kk)
      tigr_res=dble(reshape(t_rupt,resh))
      do j=1,NW(kk)
        do i=1,NL(kk)
          ruptime(i,j,kk)=interp2D(tigrNL(kk),tigrL(1:tigrNL(kk),kk),tigrNW(kk),&
           tigrW(1:tigrNW(kk),kk),tigr_res,cellL(i,kk),cellW(j,kk))&
            + m_hypo(3,kk,im)
        enddo
      enddo
      deallocate(hs,t_rupt)
      deallocate(tigr_res)
      ! In case of error
      if(iost.ne.0)then
        error=1
        return
      endif
    enddo
    
    !----------------------------------
    ! Prepare shifted rake control points (-179 is close to +179)
    allocate(m_tmp_rakes(maxval(NPL),maxval(NPW),NSeg))
    allocate(rake0(NSeg))
    do kk=1,NSeg
      rake0(kk) = m_rakes(1,1,kk,im)
      do j=1,NPW(kk)
        do i=1,NPL(kk)
          m_tmp_rakes(i,j,kk) = m_rakes(i,j,kk,im) - rake0(kk)
          if(m_tmp_rakes(i,j,kk)<-180.d0)then
            m_tmp_rakes(i,j,kk)=m_tmp_rakes(i,j,kk)+360.d0
          elseif(m_tmp_rakes(i,j,kk)>180.d0)then
            m_tmp_rakes(i,j,kk)=m_tmp_rakes(i,j,kk)-360.d0
          endif
        enddo
      enddo
    enddo
    
    !----------------------------------
    ! Interpolate control points to cell grid
    do kk=1,NSeg
      do j=1,NW(kk)
        do i=1,NL(kk)
          ! risetime
          risetime(i,j,kk)=interp2D(NPL(kk),cpL(1:NPL(kk),kk),NPW(kk),cpW(1:NPW(kk),kk),&
            m_risetime(1:NPL(kk),1:NPW(kk),kk,im),cellL(i,kk),cellW(j,kk))
          ! peaktime
          peaktime(i,j,kk)=interp2D(NPL(kk),cpL(1:NPL(kk),kk),NPW(kk),cpW(1:NPW(kk),kk),&
            m_peaktime(1:NPL(kk),1:NPW(kk),kk,im),cellL(i,kk),cellW(j,kk))
          ! rakes
          rakes(i,j,kk)=interp2D(NPL(kk),cpL(1:NPL(kk),kk),NPW(kk),cpW(1:NPW(kk),kk),&
            m_tmp_rakes(1:NPL(kk),1:NPW(kk),kk),cellL(i,kk),cellW(j,kk))
          rakes(i,j,kk) = rakes(i,j,kk) + rake0(kk)
          if(rakes(i,j,kk)<-180.d0)then
            rakes(i,j,kk)=rakes(i,j,kk)+360.d0
          elseif(rakes(i,j,kk)>180.d0)then
            rakes(i,j,kk)=rakes(i,j,kk)-360.d0
          endif
        enddo
      enddo
    enddo
    deallocate(m_tmp_rakes,rake0)
    
    ! Cycle to linearly scale total seismic moment
    iterN = 20
    do iter=1,iterN
      
      !----------------------------------
      ! Interpolate slip by bspline to cell grid
      !    (50% of the computation time of InterpMod() is here in bspline interpolation)
      do kk=1,NSeg
        tmp=splineN0(kk)+m_Nspline(kk,im)
        slip(1:NL(kk),1:NW(kk),kk)=bspline(tmp,NL(kk),NW(kk),m_slip(1:tmp,1,kk,im),&
            m_slip(1:tmp,2,kk,im),m_slip(1:tmp,3,kk,im),cellL(1:NL(kk),kk),cellW(1:NW(kk),kk))
        ! In case of error
        if(slip(1,1,kk).eq.-12345.d0)then
          error=2
          return
        endif
      enddo
      ! Remove negative slip
      where(slip.LT.0.d0) slip=0.d0
    
      !----------------------------------
      ! Generates slip-rate time function
      timefunc=0.d0
      do kk=1,NSeg
        do j=1,NW(kk)
          do i=1,NL(kk)
            ! Regularized Yoffe slip-rate time function
            timefunc(1:Ssvd,i,j,kk)=RYoffe(Ssvd,dt,ruptime(i,j,kk),risetime(i,j,kk),peaktime(i,j,kk),slip(i,j,kk))
            ! In case of error
            if(timefunc(1,i,j,kk).eq.-12345.d0)then
              error=3
              return
            endif
          enddo
        enddo
      enddo
    
      !----------------------------------
      !  Compute total scalar seismic moment
      moment=0.d0
      do kk=1,NSeg
        do k=1,Ssvd
          moment = moment + sum(timefunc(k,1:NL(kk),1:NW(kk),kk)*mu(1:NL(kk),1:NW(kk),kk))*elem(kk)*dt
        enddo
      enddo
      
      !----------------------------------
      !  Justify m_slip to reache target total seismic moment
      if(targetM0 <= 0.d0)then
        !  Disable target moment
        exit
      elseif(abs(targetM0-moment) <= 1.d9)then
        !  Seismic moment difference less then Mw0 seismic moment (success)
        exit
      elseif(iter.eq.iterN)then
        !  Target seismic moment didnt reached
        error=4
      else
        do kk=1,NSeg
          tmp=splineN0(kk)+m_Nspline(kk,im)
          m_slip(splineN0(kk)+1:tmp,3,kk,im)=m_slip(splineN0(kk)+1:tmp,3,kk,im) * (targetM0/moment)
        enddo
        cycle
      endif
      
    enddo
    
    !----------------------------------
    !  Regroup timefunc into tsou and apply rake coef. 
    coun=0
    do r=1,nrake
      do kk=1,NSeg
        do j=1,NW(kk)
          do i=1,NL(kk)
            coun=coun+1
            rangd=abs(rakes(i,j,kk)-lrake(r,kk))
            if(rangd>180.d0)then
              rangd=360.d0-rangd
            endif
            tfin(((coun-1)*Ssvd)+1:coun*Ssvd)=timefunc(1:Ssvd,i,j,kk)*dcos(rangd/180.d0*PI)
          enddo
        enddo
      enddo
    enddo
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE SaveSlipMod(fileprefix)
!---------------------------------------------------------------------
!  Save slip model into file 
!  (i.e. slip,risetime,ruptime,peaktime,rakes,tsou)
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    implicit none
    character(len=80):: fileprefix
    character(len=99):: filename
    real(8),allocatable,dimension(:):: posL,posW
    integer:: i,j,kk
    
    !  Headers of files with results
    filename=trim(fileprefix)//'_headers.txt'
    open(230,FILE=filename)
    write(230,'(A)') '#Number of segments'
    write(230,'(I5)') NSeg
    do kk=1,NSeg
      write(230,'(A8,I5)') '#Segment',kk
      allocate(posL(NL(kk)+1),posW(NW(kk)+1))
      do i=1,NL(kk)+1
        posL(i) = (dble(i-1))*leng(kk)/dble(NL(kk))
      enddo
      do j=1,NW(kk)+1
        posW(j) = (dble(j-1))*widt(kk)/dble(NW(kk))
      enddo
      write(230,'(I5)') NL(kk)+1
      write(230,'(1000F10.1)') posL
      write(230,'(I5)') NW(kk)+1
      write(230,'(1000F10.1)') posW
      deallocate(posL,posW)
    enddo
    close(230)
    
    !  1D slip-rate model as the sum of slip-rate functions over fault dip
    !filename=trim(fileprefix)//'_tfin_1D.txt'
    !open(201,FILE=filename)
    !do kk=1,NSeg
    !  do i=1,Ssvd
    !    write(201,'(1000E13.5)')(sum(timefunc(i,j,1:NW(kk),kk))*dW(kk),j=1,NL(kk))
    !  enddo
    !  write(201,*)
    !  write(201,*)
    !enddo
    !close(201)
    
    !  Full slip-rate model as one vector
    filename=trim(fileprefix)//'_tfin.txt'
    open(231,FILE=filename)
    write(231,'(1E24.15E3)') tfin
    close(231)
    
    !  Total Slip
    filename=trim(fileprefix)//'_slip2D.txt'
    open(232,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(232,'(1000E16.7E3)')(sum(timefunc(:,i,j,kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),j,kk))*dt
!       write(232,'(1000E16.7E3)')(sum(timefunc(:,i,j,kk))*dt*mu(i,j,kk)*elem(kk),i=1,NL(kk)),sum(timefunc(:,NL(kk),j,kk))*dt*mu(NL(kk),j,kk)*elem(kk)     !Output is scalar moment instead of slip
      enddo
      write(232,'(1000E16.7E3)')(sum(timefunc(:,i,NW(kk),kk))*dt,i=1,NL(kk)),sum(timefunc(:,NL(kk),NW(kk),kk))*dt
!      write(232,'(1000E16.7E3)')(sum(timefunc(:,i,NW(kk),kk))*dt*mu(i,NW,kk)*elem(kk),i=1,NL(kk)),sum(timefunc(:,NL(kk),NW(kk),kk))*dt*mu(NL(kk),NW(kk),kk)*elem(kk)   !Output is scalar moment instead of slip
      write(232,*)
      write(232,*)
    enddo
    close(232)
    
    !  Rupture time
    filename=trim(fileprefix)//'_ruptime2D.txt'
    open(233,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(233,'(1000E16.7E3)')ruptime(1:NL(kk),j,kk),ruptime(NL(kk),j,kk)
      enddo
      write(233,'(1000E16.7E3)')ruptime(1:NL(kk),NW(kk),kk),ruptime(NL(kk),NW(kk),kk)
      write(233,*)
      write(233,*)
    enddo
    close(233)
    
    !  Risetime time
    filename=trim(fileprefix)//'_risetime2D.txt'
    open(234,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(234,'(1000F8.3)')risetime(1:NL(kk),j,kk),risetime(NL(kk),j,kk)
      enddo
      write(234,'(1000F8.3)')risetime(1:NL(kk),NW(kk),kk),risetime(NL(kk),NW(kk),kk)
      write(234,*)
      write(234,*)
    enddo
    close(234)
    
    !  Peaktime time
    filename=trim(fileprefix)//'_peaktime2D.txt'
    open(235,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(235,'(1000F8.3)')peaktime(1:NL(kk),j,kk),peaktime(NL(kk),j,kk)
      enddo
      write(235,'(1000F8.3)')peaktime(1:NL(kk),NW(kk),kk),peaktime(NL(kk),NW(kk),kk)
      write(235,*)
      write(235,*)
    enddo
    close(235)
    
    !  Rakes
    filename=trim(fileprefix)//'_rakes2D.txt'
    open(236,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(236,'(1000F7.1)')rakes(1:NL(kk),j,kk),rakes(NL(kk),j,kk)
      enddo
      write(236,'(1000F7.1)')rakes(1:NL(kk),NW(kk),kk),rakes(NL(kk),NW(kk),kk)
      write(236,*)
      write(236,*)
    enddo
    close(236)
    
    !  Peak Slip-rate
    filename=trim(fileprefix)//'_sliprate2D.txt'
    open(237,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(237,'(1000E16.7E3)')(maxval(timefunc(:,i,j,kk)),i=1,NL(kk)),maxval(timefunc(:,NL(kk),j,kk))
      enddo
      write(237,'(1000E16.7E3)')(maxval(timefunc(:,i,NW(kk),kk)),i=1,NL(kk)),maxval(timefunc(:,NL(kk),NW(kk),kk))
      write(237,*)
      write(237,*)
    enddo
    close(237)
    
    !  Peak Slip-rate time
    filename=trim(fileprefix)//'_ratetime2D.txt'
    open(238,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)
        write(238,'(1000E16.7E3)')ruptime(1:NL(kk),j,kk)+peaktime(1:NL(kk),j,kk),ruptime(NL(kk),j,kk)+peaktime(NL(kk),j,kk)
      enddo
      write(238,'(1000E16.7E3)')ruptime(1:NL(kk),NW(kk),kk)+peaktime(1:NL(kk),NW(kk),kk), &
                                ruptime(NL(kk),NW(kk),kk)+peaktime(NL(kk),NW(kk),kk)
      write(238,*)
      write(238,*)
    enddo
    close(238)
    
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE Setuptempladder(nchains,tlow,thigh,dir,Chaintemp)
!---------------------------------------------------------------------
!     Setuptempladder - User routine to set up temperatures for each chain
!     Input:
!        nchains    Integer   : Number of chains
!        Tlow       Double    : Lowest temperature
!        Thigh      Double    : Highest temperature
!        dir        Character(len=80)
!     Output:
!        Chaintemp  Double array (nchains)    : Temperatures of each chain
!                                               (on this processor)
!     Notes:
!           This utility routine calculates and returns the array Chaintemps.
!           It calculates a random log-uniform set of temperature values
!           between input bounds and put the results in Chaintemps.
!           Chaintemps is written out to file `tlevels'.
!---------------------------------------------------------------------
    use ExeBox
    
#if defined MPI
    include "mpif.h"
    integer,dimension(MPI_STATUS_SIZE):: status
#endif
    integer,intent(in):: nchains
    real(8),intent(in):: tlow,thigh
    character(len=80),intent(in):: dir
    real(8):: Chaintemp(nchains)
    real(8),allocatable:: AllTemps(:)
    real(8):: aval,bval
    integer:: it,k,ierror
    Character(len=90) filename
    
    !----------------------------------
    ! Selected temperatures randomly using log-uniform distribution
    aval = log(tlow)
    bval = log(thigh)
    do it=1,nchains
      Chaintemp(it) = exp(aval + ran3(iseed)*(bval-aval))
    enddo
    
    !----------------------------------
    ! Force first chain to be at tlow
    !if(rank==1)Chaintemp(1) = tlow       ! Force some chains to be at tlow
    Chaintemp(1) = tlow
    
    !----------------------------------
    ! Write to file
#if defined MPI
    allocate(AllTemps(nchains*nproc))
    AllTemps = 0.d0
    ! Send all Temperatures to master for output to file (for diagnostics)
    call MPI_GATHER(Chaintemp,nchains,MPI_DOUBLE_PRECISION,&
        AllTemps,nchains,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    if(rank == 0)then
      filename = trim(dir)//'tlevels'
      open(15,file=filename,status='unknown')
      k = nchains*(nproc-1)
      if(k.gt.1)call HPSORT(k,AllTemps(nchains+1)) ! Order the temperatures for neat output
      do it = nchains+1,nproc*nchains
        write(15,*)it-nchains,AllTemps(it)
      enddo
      close(15)
    endif
    deallocate(AllTemps)
#else
    filename = trim(dir)//'log/tlevels.log'
    open(15,file=filename,status='unknown')
    do it=1,nchains
      write(15,*)it,Chaintemp(it)
    enddo
    close(15)
#endif
    RETURN
!---------------------------------------------------------------------
END SUBROUTINE


