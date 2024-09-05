!---------------------------------------------------------------------
!
!  Execute PSI after processing
!
!---------------------------------------------------------------------
!
!  Author: Miroslav Hallo
!  Charles University, Faculty of Mathematics and Physics
!  E-mail: hallo@karel.troja.mff.cuni.cz
!  Revision 10/2017: The initial version
!  Revision 2018: Fixed multiple minor issues
!  Revision 2019: Added REAL4
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

MODULE FinBox
!---------------------------------------------------------------------
!  Module containing global after processing veriables
!---------------------------------------------------------------------
    ! Global slip model bins
    INTEGER:: nBins
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: allBins
    ! Global slip model histograms
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:):: pop_slip,pop_risetime,pop_ruptime
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:):: pop_peaktime,pop_rakes
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:):: pop_sliprate,pop_ratetime
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:):: pop_joint,cou_hypo
    ! Rupture velocity
    REAL(8),ALLOCATABLE,DIMENSION(:,:,:):: ruptveloc
    ! Vars for statistics
    REAL(8):: MaxSlip,MaxTime,MaxMisf,MaxSlipRate,wam(2,4),MinVR
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: cellLm,cellWm
    REAL(8),ALLOCATABLE,DIMENSION(:):: MeanRupVel,BestRupVel
    INTEGER:: cellmult
    ! The best model
    REAL(8):: best_misfit,map_misfit
    ! Model count
    INTEGER:: mcount
    ! MAP from N splines
    INTEGER,ALLOCATABLE,DIMENSION(:):: map_Nsplines
    ! NeoBins for number of spline points
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: neoBins
    INTEGER,ALLOCATABLE,DIMENSION(:,:):: pop_neo
    ! SPECIAL FOR SLIP-RATE SLICES
    INTEGER:: srsample,srs
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:):: pop_sliprate1,pop_sliprate2,pop_sliprate3,pop_sliprate4,pop_sliprate5
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:):: pop_sliprate6,pop_sliprate7,pop_sliprate8,pop_sliprate9
!---------------------------------------------------------------------
END MODULE



SUBROUTINE FinExe()
!---------------------------------------------------------------------
!  Creates initial model parameters and execute inversion,
!  the model space is prospected by a Monte-Carlo method
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    use LibBox
    use FinBox
    implicit none
    integer:: run,i,j,kk,p,ios,popInd,error,im,imB,imM,imOUT,ncount
    character(len=80):: dir
    character(len=96):: filename
    character(len=80):: fileprefix
    character(len=4):: cproc,mfnam
    real(8):: WL
#ifdef REAL4
    real(4),allocatable:: finUvec(:),Hnew(:,:)
#else
    real(8),allocatable:: finUvec(:),Hnew(:,:)
#endif
    real(8):: logPPD,misfit
    logical :: file_exists
    
    !----------------------------------
    !  Read parameters from afterprocess file
	open(17,file='input/afterprocess.dat',action='read',status='old')
    read(17,*)
    read(17,*) nBins
    read(17,*)
    read(17,*) MaxSlip
    read(17,*)
    read(17,*) MaxTime
    read(17,*)
    read(17,*) MaxSlipRate
    read(17,*)
    read(17,*) MinVR
    read(17,*)
    read(17,*) cellmult
    read(17,*)
    read(17,*) srsample
    close(17)
	
    dir = './inv/' ! work directory with models and for outputs
    mfnam='00' ! text identificator of the output files
    nproc = 1000 ! maximal number of CPU procesors (if you add more, do not worry non-existing files are skiped)
	
    !nBins = 23 ! number of bins for statistics
    !MaxSlip = 8.0d0 ! Maximum expected slip [m] (for bins)
    !MaxTime = 20.d0 ! Maximum expected rupture time [s] (for bins)
    !MaxSlipRate = 5.0d0 ! Maximum expected peak slip-rate [m/s] (for bins)
    !MinVR = 70.d0 ! Minimum data variance reduction (for MaxMisf and misfit bins)  
    !cellmult = 5 ! oversampling cell grid number for hypocentrum stats
    !srsample = 5 ! Each Xth slip-rate function sample for creation of the time slice (*dt for [s])
    
    !----------------------------------
    ! Initial arrays for after processing
    write(*,*) ' Initial arrays'
    call InitFin()
    allocate(map_Nsplines(NSeg))
    
    !----------------------------------
    ! Loop over search run
    write(*,*) ' Process all models'
    im=1 ! test model index
    imB=2 ! The best model index (max likelyhood from all)
    imM=3 ! The MAP model index (max likelyhood from MAP dimension)
    
    do run=1,2
      ! Zeros into allocated population arrays
      pop_slip=0
      pop_risetime=0
      pop_ruptime=0
      pop_peaktime=0
      pop_rakes=0
      pop_joint=0
      cou_hypo=0
      pop_sliprate=0
      pop_ratetime=0
      pop_neo=0
      
      pop_sliprate1=0
      pop_sliprate2=0
      pop_sliprate3=0
      pop_sliprate4=0
      pop_sliprate5=0
      pop_sliprate6=0
      pop_sliprate7=0
      pop_sliprate8=0
      pop_sliprate9=0
      
      best_misfit = 1d10 ! initial (very high) misfit
      map_misfit = 1d10 ! initial (very high) misfit
      mcount=0
      
      if(run.eq.1)then
        write(*,*) ' THE FIRST RUN - POP and BEST'
      else
        write(*,*) ' THE SECOND RUN - MAP'
      endif
      
      !----------------------------------
      ! Loop over all files with results
      do p=0,nproc-1
        ! File
        write(cproc,'(I0)') p
        filename=trim(dir)//'xmodels_'//trim(cproc)//'.dat'
        ! Check if the file exists
        inquire(file=filename, exist=file_exists)
        if(.NOT.file_exists)then
          cycle
        endif
        
        ! Open the file
        ifile = 73
        open(unit=ifile,form='unformatted',file=filename,action='read')
        write(*,'(A98)') filename
      
        !----------------------------------
        ! Loop over models within the file
        do while(.true.)
          ! Read new model
          read(ifile,iostat=ios) logPPD,misfit
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_veloc(:,:,:,im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_risetime(:,:,:,im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_peaktime(:,:,:,im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_rakes(:,:,:,im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_hypo(:,:,im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_Nspline(:,im)
          if(ios.ne.0) exit 
          read(ifile,iostat=ios) m_slip(:,:,:,im)
          if(ios.ne.0) exit
          
          mcount=mcount+1
          
          ! Screen output
          if(mod(mcount,1000).eq.0)then
            write(*,'(I11,A17)') mcount,' models processed'
          endif
          
          ! Model parameters into slip,risetime,ruptime,peaktime,rakes
          call InterpMod(im,-12345.d0,error)
          if(error.ne.0)then
            write(*,*) ' Error in the interpretation of model parameters!',error
            stop
            !cycle
          endif
          
          ! Prevent zero division in weighted averages
          if(moment.le.0.d0)then
            do kk=1,NSeg
              do j=1,NW(kk)
                slip(1:NL(kk),j,kk) = 1.d-99
              enddo
            enddo
          endif
          
          ! Interpolate rupture velocity control points to cell grid
          !----------------------------------
          do kk=1,NSeg
            do j=1,NW(kk)
              do i=1,NL(kk)
              ! rupture velocity
              ruptveloc(i,j,kk)=interp2D(NPL(kk),cpL(1:NPL(kk),kk),NPW(kk),cpW(1:NPW(kk),kk),&
                m_veloc(1:NPL(kk),1:NPW(kk),kk,im),cellL(i,kk),cellW(j,kk))
              enddo
            enddo
          enddo
          
          ! Count population
          if(run.eq.1)then
            do kk=1,NSeg
              wam = 0.d0
              do j=1,NW(kk)
                do i=1,NL(kk)
                  ! Slip
                  popInd = binarysearch(nBins,allBins(:,1,kk),slip(i,j,kk))
                  pop_slip(i,j,kk,popInd) = pop_slip(i,j,kk,popInd)+1
                  ! Risetime
                  popInd = binarysearch(nBins,allBins(:,2,kk),risetime(i,j,kk))
                  pop_risetime(i,j,kk,popInd) = pop_risetime(i,j,kk,popInd)+1
                  wam(1,1) = wam(1,1) + risetime(i,j,kk)*slip(i,j,kk) ! Weighted mean
                  wam(2,1) = wam(2,1) + slip(i,j,kk) ! Weight
                  ! Ruptime
                  popInd = binarysearch(nBins,allBins(:,3,kk),ruptime(i,j,kk))
                  pop_ruptime(i,j,kk,popInd) = pop_ruptime(i,j,kk,popInd)+1
                  ! Peaktime
                  popInd = binarysearch(nBins,allBins(:,4,kk),peaktime(i,j,kk))
                  pop_peaktime(i,j,kk,popInd) = pop_peaktime(i,j,kk,popInd)+1
                  wam(1,2) = wam(1,2) + peaktime(i,j,kk)*slip(i,j,kk) ! Weighted mean
                  wam(2,2) = wam(2,2) + slip(i,j,kk) ! Weight
                  ! Rakes
                  if(rake(kk)>90.d0 .and. rakes(i,j,kk)<-90.d0)then
                    popInd = binarysearch(nBins,allBins(:,5,kk),(rakes(i,j,kk)+360.d0))
                    wam(1,3) = wam(1,3) + (rakes(i,j,kk)+360.d0)*slip(i,j,kk) ! Weighted mean
                    wam(2,3) = wam(2,3) + slip(i,j,kk) ! Weight
                  elseif(rake(kk)<-90.d0 .and. rakes(i,j,kk)>90.d0)then
                    popInd = binarysearch(nBins,allBins(:,5,kk),(rakes(i,j,kk)-360.d0))
                    wam(1,3) = wam(1,3) + (rakes(i,j,kk)-360.d0)*slip(i,j,kk) ! Weighted mean
                    wam(2,3) = wam(2,3) + slip(i,j,kk) ! Weight
                  else
                    popInd = binarysearch(nBins,allBins(:,5,kk),rakes(i,j,kk))
                    wam(1,3) = wam(1,3) + rakes(i,j,kk)*slip(i,j,kk) ! Weighted mean
                    wam(2,3) = wam(2,3) + slip(i,j,kk) ! Weight
                  endif
                  pop_rakes(i,j,kk,popInd) = pop_rakes(i,j,kk,popInd)+1
                  ! Peak Slip-Rate
                  popInd = binarysearch(nBins,allBins(:,11,kk),maxval(timefunc(:,i,j,kk)))
                  pop_sliprate(i,j,kk,popInd) = pop_sliprate(i,j,kk,popInd)+1
                  ! Peak Slip-Rate slices
				  srs = srsample+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate1(i,j,kk,popInd) = pop_sliprate1(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*2)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate2(i,j,kk,popInd) = pop_sliprate2(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*3)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate3(i,j,kk,popInd) = pop_sliprate3(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*4)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate4(i,j,kk,popInd) = pop_sliprate4(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*5)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate5(i,j,kk,popInd) = pop_sliprate5(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*6)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate6(i,j,kk,popInd) = pop_sliprate6(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*7)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate7(i,j,kk,popInd) = pop_sliprate7(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*8)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate8(i,j,kk,popInd) = pop_sliprate8(i,j,kk,popInd)+1
                  endif
				  srs = (srsample*9)+1
                  if(srs.le.Ssvd)then
                    popInd = binarysearch(nBins,allBins(:,11,kk),timefunc(srs,i,j,kk))
                    pop_sliprate9(i,j,kk,popInd) = pop_sliprate9(i,j,kk,popInd)+1
                  endif
                  
                  ! Peak Slip-Rate time
                  popInd = binarysearch(nBins,allBins(:,12,kk),ruptime(i,j,kk)+peaktime(i,j,kk))
                  pop_ratetime(i,j,kk,popInd) = pop_ratetime(i,j,kk,popInd)+1
                  ! Rupture velocity
                  wam(1,4) = wam(1,4) + ruptveloc(i,j,kk)*slip(i,j,kk) ! Weighted mean
                  wam(2,4) = wam(2,4) + slip(i,j,kk) ! Weight
                enddo
              enddo
              ! Segment delay 
              popInd = binarysearch(nBins,allBins(:,6,kk),m_hypo(3,kk,im))
              pop_joint(1,kk,popInd) = pop_joint(1,kk,popInd)+1
              ! N splines
              popInd = binarysearch(nBins,allBins(:,7,kk),dble(m_Nspline(kk,im)))
              pop_joint(2,kk,popInd) = pop_joint(2,kk,popInd)+1
              ! Misfit
              popInd = binarysearch(nBins,allBins(:,8,kk),misfit)
              pop_joint(3,kk,popInd) = pop_joint(3,kk,popInd)+1
              ! Total moment
              popInd = binarysearch(nBins,allBins(:,9,kk),moment)
              pop_joint(4,kk,popInd) = pop_joint(4,kk,popInd)+1
              ! Magnitude
              popInd = binarysearch(nBins,allBins(:,10,kk),((log10(max(1.d0,moment))-9.1d0)/1.5d0))
              pop_joint(5,kk,popInd) = pop_joint(5,kk,popInd)+1
              ! Risetime (Weighted arithmetic mean)
              popInd = binarysearch(nBins,allBins(:,2,kk),(wam(1,1)/wam(2,1)))
              pop_joint(6,kk,popInd) = pop_joint(6,kk,popInd)+1
              ! Peaktime (Weighted arithmetic mean)
              popInd = binarysearch(nBins,allBins(:,4,kk),(wam(1,2)/wam(2,2)))
              pop_joint(7,kk,popInd) = pop_joint(7,kk,popInd)+1
              ! Rakes (Weighted arithmetic mean)
              popInd = binarysearch(nBins,allBins(:,5,kk),(wam(1,3)/wam(2,3)))
              pop_joint(8,kk,popInd) = pop_joint(8,kk,popInd)+1
              ! Rupture velocity (Weighted arithmetic mean)
              MeanRupVel(kk) = wam(1,4)/wam(2,4)
              popInd = binarysearch(nBins,allBins(:,13,kk),(wam(1,4)/wam(2,4)))
              pop_joint(9,kk,popInd) = pop_joint(9,kk,popInd)+1
              
              ! Hypocenter position
              i = binarysearch(NL(kk)*cellmult,cellLm(1:NL(kk)*cellmult,kk),m_hypo(1,kk,im))
              j = binarysearch(NW(kk)*cellmult,cellWm(1:NW(kk)*cellmult,kk),m_hypo(2,kk,im))
              cou_hypo(i,j,kk) = cou_hypo(i,j,kk)+1
              
              ! NeoBins for N splines
              popInd = binarysearch(MaxNSpline,neoBins(:,kk),dble(m_Nspline(kk,im)))
              pop_neo(kk,popInd) = pop_neo(kk,popInd)+1
            enddo
            
          endif
          
          ! The best model or MAP
          if(run.eq.1)then   ! The best model
            if(misfit.lt.best_misfit)then
              m_veloc(:,:,:,imB) = m_veloc(:,:,:,im)
              m_risetime(:,:,:,imB) = m_risetime(:,:,:,im)
              m_peaktime(:,:,:,imB) = m_peaktime(:,:,:,im)
              m_rakes(:,:,:,imB) = m_rakes(:,:,:,im)
              m_hypo(:,:,imB) = m_hypo(:,:,im)
              m_Nspline(:,imB) = m_Nspline(:,im)
              m_slip(:,:,:,imB) = m_slip(:,:,:,im)
              best_misfit = misfit
              BestRupVel(:) = MeanRupVel(:)
              !write(*,*) misfit
            endif
          else   ! The MAP model
            ncount = 0
            ! Find N splines of current model
            do kk=1,NSeg
              !popInd = binarysearch(nBins,allBins(:,7,kk),dble(m_Nspline(kk,im)))
              popInd = binarysearch(MaxNSpline,neoBins(:,kk),dble(m_Nspline(kk,im)))
              if(popInd.eq.map_Nsplines(kk))then
                ncount = ncount+1
              endif
            enddo
            if(ncount.eq.NSeg)then
              if(misfit.lt.map_misfit)then
                m_veloc(:,:,:,imM) = m_veloc(:,:,:,im)
                m_risetime(:,:,:,imM) = m_risetime(:,:,:,im)
                m_peaktime(:,:,:,imM) = m_peaktime(:,:,:,im)
                m_rakes(:,:,:,imM) = m_rakes(:,:,:,im)
                m_hypo(:,:,imM) = m_hypo(:,:,im)
                m_Nspline(:,imM) = m_Nspline(:,im)
                m_slip(:,:,:,imM) = m_slip(:,:,:,im)
                map_misfit = misfit
                BestRupVel(:) = MeanRupVel(:)
                !write(*,*) misfit
              endif
            endif
          endif
        
        enddo  ! End loop over models within the file
      
        !----------------------------------
        ! Close the file
        close(ifile)
      enddo ! End loop over all files
    
      !----------------------------------
      ! Prepare variables after the first run
      if(run.eq.1)then
        ! Find spline N of MAP for all segments
        do kk=1,NSeg
          !map_Nsplines(kk) = maxloc(pop_joint(2,kk,1:nBins),1)
          map_Nsplines(kk) = maxloc(pop_neo(kk,1:MaxNSpline),1)
        enddo
        imOUT = imB
      else
        imOUT = imM
      endif
      
      !----------------------------------
      ! Write stats of after processing
      if(run.eq.1)then
        fileprefix = 'inv/pop'//trim(mfnam)
        call SaveStats(fileprefix)
      endif
      
      !----------------------------------
      ! Interpret and write the best model
      call Forward(imOUT,-12345.d0,misfit,logPPD,error)
      if(run.eq.1)then
        fileprefix = 'inv/best'//trim(mfnam)
      else
        fileprefix = 'inv/map'//trim(mfnam)
      endif
      call SaveSlipMod(fileprefix)
      
      !----------------------------------
      ! Save final the best fitting synthetic waveforms
      allocate(finUvec(Nseis))
#ifdef MKL
#ifdef REAL4
      call sgemv('N',Nseis,Msvd,1.,H,Nseis,tfin,1,0.,finUvec,1)
#else 
      call dgemv('N',Nseis,Msvd,1.d0,H,Nseis,tfin,1,0.d0,finUvec,1)
#endif
#else
      finUvec=matmul(H,tfin)
#endif
      
      ! Write standardized waveforms
      filename=trim(fileprefix)//'_Unezstd.txt'
      open(232,FILE=filename)
      do i=1,nT
        write(232,'(1000E17.7E3)')dt*(iT1-1+i-1),(finUvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(232)
      
      ! Read original H from file and create non-standardized waveforms
      allocate(Hnew(Nseis,Msvd))
      open(422,FILE='var_H.tmp',ACCESS='stream',FORM='unformatted')
      do i=1,Nseis
        read(422) Hnew(i,1:Msvd)
      enddo
      close(422)
#ifdef MKL
#ifdef REAL4
      call sgemv('N',Nseis,Msvd,1.,Hnew,Nseis,tfin,1,0.,finUvec,1)
#else
      call dgemv('N',Nseis,Msvd,1.d0,Hnew,Nseis,tfin,1,0.d0,finUvec,1)
#endif
#else
      finUvec=matmul(Hnew,tfin)
#endif
      
      ! Write non-standardized waveforms
      filename=trim(fileprefix)//'_Unez.txt'
      open(239,FILE=filename)
      do i=1,nT
        write(239,'(1000E17.7E3)')dt*(iT1-1+i-1),(finUvec((j-1)*nT+i),j=1,NSTAcomp)
      enddo
      close(239)
      
      deallocate(Hnew)
      deallocate(finUvec)
    
      !----------------------------------
      !  Write joint the best model features
      filename=trim(fileprefix)//'_joint.txt'
      open(236,FILE=filename)
      write(236,'(A)') '#Hypocenter L,W,delay [m,m,s]'
      write(236,'(100F10.1)') m_hypo(1,1:NSeg,imOUT)
      write(236,'(100F10.1)') m_hypo(2,1:NSeg,imOUT)
      write(236,'(100F10.4)') m_hypo(3,1:NSeg,imOUT)
      write(236,'(A)') '#Spline points number'
      write(236,'(100I10)') m_Nspline(1:NSeg,imOUT)
      write(236,'(A)') '#Average rupture spead'
      write(236,'(100F10.1)') BestRupVel(1:NSeg)
      write(236,'(A)') '#Total scalar moment'
      write(236,'(E10.3)') moment
      write(236,'(A)') '#Moment magnitude'
      write(236,'(F10.3)') (log10(moment)-9.1d0)/1.5d0
      write(236,'(A)') '#Misfit'
      write(236,'(F10.3)') misfit
      write(236,'(A)') '#Variance reduction'
      write(236,'(F10.3)')  (1-(misfit/sqUvec))*100.d0
      write(236,'(A)') '#Number of similar models'
      write(236,'(I11)') mcount
      close(236)
      
      !  Save the best model parameters
      filename=trim(fileprefix)//'_model.txt'
      open(237,FILE=filename)
      do kk=1,NSeg
        write(237,'(A,I5)') '#Segment',kk
        write(237,'(A)') '#Velocity [m/s]'
        do j=1,NPW(kk)
          write(237,'(1000F8.1)') m_veloc(1:NPL(kk),j,kk,imOUT)
        enddo
        write(237,'(A)') '#Risetime [s]'
        do j=1,NPW(kk)
          write(237,'(1000F8.2)') m_risetime(1:NPL(kk),j,kk,imOUT)
        enddo
        write(237,'(A)') '#Peaktime [s]'
        do j=1,NPW(kk)
          write(237,'(1000F8.3)') m_peaktime(1:NPL(kk),j,kk,imOUT)
        enddo
        write(237,'(A)') '#Rakes [deg]'
        do j=1,NPW(kk)
          write(237,'(1000F8.1)') m_rakes(1:NPL(kk),j,kk,imOUT)
        enddo
        write(237,'(A)') '#Hypocenter L,W,delay [m,m,s]'
        write(237,'(F9.1,F9.1,F9.3)') m_hypo(1,kk,imOUT), m_hypo(2,kk,imOUT), m_hypo(3,kk,imOUT)
        write(237,'(A)') '#Spline points number (total,active)'
        write(237,'(2I10)') m_Nspline(kk,imOUT)+splineN0(kk), m_Nspline(kk,imOUT)
        write(237,'(A)') '#Spline points L,W,val [m,m,m]'
        write(237,'(1000F9.1)') m_slip(1:m_Nspline(kk,imOUT)+splineN0(kk),1,kk,imOUT)
        write(237,'(1000F9.1)') m_slip(1:m_Nspline(kk,imOUT)+splineN0(kk),2,kk,imOUT)
        write(237,'(1000F9.3)') m_slip(1:m_Nspline(kk,imOUT)+splineN0(kk),3,kk,imOUT)
      enddo
      close(237)
      
    enddo
    
!---------------------------------------------------------------------
CONTAINS



SUBROUTINE InitFin()
!---------------------------------------------------------------------
!  Initial after processing
!---------------------------------------------------------------------
    implicit none
    integer kk,b,im,pert
    
    pert = 3    ! 1=test, 2=best, 3=MAP
    sqUvec = dot_product(Uvec,Uvec)
    
    !----------------------------------
    ! Maximum expected misfit (for bins)
    MaxMisf = sqUvec*(1.d0 - MinVR/100.d0)
    
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
    allocate(m_veloc(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_risetime(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_peaktime(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_rakes(maxval(NPL),maxval(NPW),NSeg,pert))
    allocate(m_hypo(3,NSeg,pert))
    allocate(m_slip(maxval(splineN0)+MaxNSpline,3,NSeg,pert))
    allocate(m_Nspline(NSeg,pert))
    
    !----------------------------------
    ! Spline zero control points on edges
    m_slip=0.d0 ! Required for zero control points on edges
    do im=1,pert ! Set for all including perturbed model
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
    ! Alocate arrays for after processing
    allocate(allBins(nBins,13,NSeg))
    allocate(pop_slip(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_risetime(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_ruptime(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_peaktime(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_rakes(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_joint(9,NSeg,nBins))
    allocate(cou_hypo(maxval(NL)*cellmult,maxval(NW)*cellmult,NSeg))
    allocate(pop_sliprate(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_ratetime(maxval(NL),maxval(NW),NSeg,nBins))
    
    ! SPECIAL slip-rate pop
    allocate(pop_sliprate1(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate2(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate3(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate4(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate5(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate6(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate7(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate8(maxval(NL),maxval(NW),NSeg,nBins))
    allocate(pop_sliprate9(maxval(NL),maxval(NW),NSeg,nBins))
	
    !----------------------------------
    ! Prepare bins for parameters
    allBins=0.d0
    do b=1,nBins
      do kk=1,NSeg
        allBins(b,1,kk) = (b-1)*MaxSlip/nBins ! Slip
        allBins(b,2,kk) = thr_rise(1) + (b-1)*(thr_rise(2)-thr_rise(1))/nBins ! Risetime
        allBins(b,3,kk) = (b-1)*MaxTime/nBins ! Ruptime
        allBins(b,4,kk) = thr_peak(1) + (b-1)*(thr_peak(2)-thr_peak(1))/nBins ! Peaktime
        allBins(b,5,kk) = rake(kk)-thr_rake(1) + (b-1)*(thr_rake(1)*2.d0)/nBins ! Rakes
        ! pop_joint
        if(kk.eq.1)then
          allBins(b,6,kk) = segDelay(kk)-thr_delay(1) + (b-1)*(thr_delay(1)*2.d0)/nBins ! 1st Segment delay
        else
          allBins(b,6,kk) = segDelay(kk)-thr_delay(2) + (b-1)*(thr_delay(2)*2.d0)/nBins ! Segment delay
        endif
        allBins(b,7,kk) = dble(b-1) ! N splines
        allBins(b,8,kk) = (b-1)*MaxMisf/nBins ! Misfit
        allBins(b,9,kk) = Mfix-thr_M0 + (b-1)*(thr_M0*2.d0)/nBins! Total moment
        allBins(b,10,kk) = (log10(max(1.d0,Mfix-thr_M0))-9.1d0)/1.5d0 + &
               (b-1)*(((log10(Mfix+thr_M0)-9.1d0)/1.5d0)-((log10(max(1.d0,Mfix-thr_M0))-9.1d0)/1.5d0))/nBins! Magnitude
        allBins(b,11,kk) = (b-1)*MaxSlipRate/nBins ! Peak Slip-Rate
        allBins(b,12,kk) = (b-1)*(MaxTime+thr_peak(2))/nBins ! Peak Slip-Rate time
        allBins(b,13,kk) = thr_velo(1) + (b-1)*(thr_velo(2)-thr_velo(1))/nBins ! Velocity
      enddo
    enddo
    
    !----------------------------------
    ! Prepare NeoBins for N splines
    allocate(neoBins(MaxNSpline,NSeg))
    allocate(pop_neo(NSeg,MaxNSpline))
    neoBins=0.d0
    do b=1,MaxNSpline
      do kk=1,NSeg
        neoBins(b,kk) = dble(b) ! N splines
      enddo
    enddo
    
    !----------------------------------
    ! Prepare oversampled cell grid
    allocate(cellLm(maxval(NL)*cellmult,NSeg),cellWm(maxval(NW)*cellmult,NSeg))
    do kk=1,NSeg
      do i=1,NL(kk)*cellmult
        cellLm(i,kk)=(dble(i-1))*dL(kk)/cellmult
      enddo
      do j=1,NW(kk)*cellmult
        cellWm(j,kk)=(dble(j-1))*dW(kk)/cellmult
      enddo
    enddo
    
    !----------------------------------
    ! Prepare rupture velocity field
    allocate(ruptveloc(maxval(NL),maxval(NW),NSeg))
    allocate(MeanRupVel(NSeg))
    allocate(BestRupVel(NSeg))
    ruptveloc=0.d0
    MeanRupVel=0.d0
    BestRupVel=0.d0
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE SaveStats(fileprefix)
!---------------------------------------------------------------------
!  Save slip model into file 
!  (i.e. slip,risetime,ruptime,peaktime,rakes,tsou)
!---------------------------------------------------------------------
    implicit none
    character(len=80):: fileprefix
    character(len=99):: filename
    real(8),allocatable,dimension(:):: posL,posW
    integer:: i,j,kk,p
    
    !  Headers of files with results
    filename=trim(fileprefix)//'_headers.txt'
    open(230,FILE=filename)
    write(230,'(A)') '#Number of segments / bins / sqUvec / models / dummy'
    write(230,'(I5)') NSeg
    write(230,'(I5)') nBins
    write(230,'(E16.7E3)') sqUvec
    write(230,'(I11)') mcount
    write(230,'(F6.3)') 0.12345d0
    do kk=1,NSeg
      write(230,'(A8,I5)') '#Segment',kk
      write(230,'(I5)') NL(kk)
      write(230,'(I5)') NW(kk)
      write(230,'(1000E16.7E3)') allBins(1:nBins,1,kk) ! Slip
      write(230,'(1000F8.3)') allBins(1:nBins,2,kk) ! Risetime
      write(230,'(1000E16.7E3)') allBins(1:nBins,3,kk) ! Ruptime
      write(230,'(1000F8.3)') allBins(1:nBins,4,kk) ! Peaktime
      write(230,'(1000F7.1)') allBins(1:nBins,5,kk) ! Rakes
      write(230,'(1000F8.3)') allBins(1:nBins,6,kk) ! Segment delay
      write(230,'(1000I10)') int(allBins(1:nBins,7,kk)) ! N splines
      write(230,'(1000E16.7E3)') allBins(1:nBins,8,kk) ! Misfit
      write(230,'(1000E16.7E3)') allBins(1:nBins,9,kk) ! Total moment
      write(230,'(1000F8.3)') allBins(1:nBins,10,kk) ! Magnitude
      write(230,'(1000E16.7E3)') allBins(1:nBins,11,kk) ! Peak Slip-Rate
      write(230,'(1000E16.7E3)') allBins(1:nBins,12,kk) ! Peak Slip-Rate time
      write(230,'(1000F7.1)') allBins(1:nBins,13,kk) ! Rupture velocity
      write(230,'(I5)') NL(kk)*cellmult
      write(230,'(I5)') NW(kk)*cellmult
      write(230,'(1000E16.7E3)') cellLm(1:NL(kk)*cellmult,kk)
      write(230,'(1000E16.7E3)') cellWm(1:NW(kk)*cellmult,kk)
      write(230,'(I5)') NL(kk)+1
      write(230,'(I5)') NW(kk)+1
      allocate(posL(NL(kk)+1),posW(NW(kk)+1))
      do i=1,NL(kk)+1
        posL(i) = (dble(i-1))*leng(kk)/dble(NL(kk))
      enddo
      do j=1,NW(kk)+1
        posW(j) = (dble(j-1))*widt(kk)/dble(NW(kk))
      enddo
      write(230,'(1000F10.1)') posL
      write(230,'(1000F10.1)') posW
      deallocate(posL,posW)
      write(230,'(I5)') NPL(kk)
      write(230,'(I5)') NPW(kk)
      write(230,'(1000F10.1)') cpL(1:NPL(kk),kk)
      write(230,'(1000F10.1)') cpW(1:NPW(kk),kk)
    enddo
    close(230)
    
    !  Total Slip
    filename=trim(fileprefix)//'_slip2D.txt'
    open(231,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(231,'(I4)') p
        do j=1,NW(kk)
          write(231,'(1000I11)') pop_slip(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(231,*)
      write(231,*)
    enddo
    close(231)
    
    !  Rupture time
    filename=trim(fileprefix)//'_ruptime2D.txt'
    open(232,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(232,'(I4)') p
        do j=1,NW(kk)
          write(232,'(1000I11)') pop_ruptime(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(232,*)
      write(232,*)
    enddo
    close(232)
    
    !  Risetime time
    filename=trim(fileprefix)//'_risetime2D.txt'
    open(233,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(233,'(I4)') p
        do j=1,NW(kk)
          write(233,'(1000I11)') pop_risetime(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(233,*)
      write(233,*)
    enddo
    close(233)
    
    !  Peaktime time
    filename=trim(fileprefix)//'_peaktime2D.txt'
    open(234,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(234,'(I4)') p
        do j=1,NW(kk)
          write(234,'(1000I11)') pop_peaktime(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(234,*)
      write(234,*)
    enddo
    close(234)
    
    !  Rakes
    filename=trim(fileprefix)//'_rakes2D.txt'
    open(235,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(235,'(I4)') p
        do j=1,NW(kk)
          write(235,'(1000I11)') pop_rakes(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(235,*)
      write(235,*)
    enddo
    close(235)
    
    !  Peak Slip-Rate
    filename=trim(fileprefix)//'_sliprate2D.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    
    !  Peak Slip-Rate slices
    filename=trim(fileprefix)//'_sliprate2D_01.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate1(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_02.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate2(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_03.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate3(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_04.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate4(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_05.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate5(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_06.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate6(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_07.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate7(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_08.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate8(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
    filename=trim(fileprefix)//'_sliprate2D_09.txt'
    open(245,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(245,'(I4)') p
        do j=1,NW(kk)
          write(245,'(1000I11)') pop_sliprate9(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(245,*)
      write(245,*)
    enddo
    close(245)
	
    !  Peak Slip-Rate time
    filename=trim(fileprefix)//'_ratetime2D.txt'
    open(255,FILE=filename)
    do kk=1,NSeg
      do p=1,nBins
        write(255,'(I4)') p
        do j=1,NW(kk)
          write(255,'(1000I11)') pop_ratetime(1:NL(kk),j,kk,p)
        enddo
      enddo
      write(255,*)
      write(255,*)
    enddo
    close(255)
    
    !  Segment delay / N splines / Misfit / Mw / Average Risetime / Average Peaktime / Average Rake / Average Velocity
    filename=trim(fileprefix)//'_joint.txt'
    open(236,FILE=filename)
    do kk=1,NSeg
      write(236,'(1000I11)') pop_joint(1,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(2,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(3,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(4,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(5,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(6,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(7,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(8,kk,1:nBins)
      write(236,'(1000I11)') pop_joint(9,kk,1:nBins)
      write(236,*)
      write(236,*)
    enddo
    close(236)
    
    !  Hypocenter position
    filename=trim(fileprefix)//'_hypo2D.txt'
    open(237,FILE=filename)
    do kk=1,NSeg
      do j=1,NW(kk)*cellmult
        write(237,'(1000I11)') cou_hypo(1:NL(kk)*cellmult,j,kk)
      enddo
      write(237,*)
      write(237,*)
    enddo
    close(237)
    
    !  NeoBins for N splines
    filename=trim(fileprefix)//'_Nsplines.txt'
    open(238,FILE=filename)
    write(238,'(A)') '#Maximum N splines / models'
    write(238,'(I5)') MaxNSpline
    write(238,'(I11)') mcount
    do kk=1,NSeg
      write(238,'(A8,I5)') '#Segment',kk
      write(238,'(1000I11)') int(neoBins(1:MaxNSpline,kk)) ! N splines
      write(238,'(1000I11)') pop_neo(kk,1:MaxNSpline)
    enddo
    close(238)
    
!---------------------------------------------------------------------
END SUBROUTINE



!---------------------------------------------------------------------
END !FinExe()

