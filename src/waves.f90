!---------------------------------------------------------------------
!
!  Read inputs, filter, integrate and cut the real waveforms for PSI (Uvec)
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
!  Revision 1/2018: The initial version
!
!  Copyright (C) 2018  Miroslav Hallo
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

PROGRAM WAVES
!---------------------------------------------------------------------
!  Main program
!  Prepare file with saved Uvec array for PSI
!---------------------------------------------------------------------
    use InputBox
    implicit none
    character(len=80):: filename
    character(len=4):: jjtxt
    real(8):: wavDT
    real(4):: dtr4,f1r4,f2r4
    real(4),allocatable,dimension(:):: fltr4,waveTmp
    real(8),allocatable:: Uvec(:)
    integer:: i,j,jj,k,mm,wavN,dspl,cou
    
    write(*,*) "------------------------------------------------------"
    write(*,*) "Read, filter, integrate and cut the real waveforms"
    !----------------------------------
    ! Read input files and init input variables
    write(*,*) "------------------------------------------------------"
    write(*,*) "Initiation of the input"
    call InitInput()
    
    write(*,*) " Number of target (synth.) time samples of one waveform:"
    write(*,'(I6)') nT
    write(*,*) " Target (synth.) dt [sec]:"
    write(*,'(F10.4)') dt
    write(*,*) " Component usage (N E Z)  Filter_freq [Hz]"
    do i=1,NRseis
      write(*,'(3I3,2F8.4)') stainfo(:,i),fc1(fcsta(i)),fc2(fcsta(i))
    enddo
    
    ! Allocate resultant data vector
    allocate(Uvec(Nseis))
    
    !----------------------------------
    ! Read waveforms ad process
    write(*,*) "------------------------------------------------------"
    write(*,*) "Read waveforms (NEZwave.dat) and process"
    open(10,file='NEZwave.dat',action='read',status='old')
    
    j=0
    do jj=1,NRseis
      write(*,*) ' Station: ',jj
      f1r4=real(fc1(fcsta(jj)))
      f2r4=real(fc2(fcsta(jj)))
      do k=1,3
        ! Skip if zero in station info
        if(stainfo(k,jj)==0)then
          write(*,*) '   -> component skipped'
          read(10,*)
          read(10,*)
          cycle
        endif
        j=j+1
        
        ! Read waveform
        read(10,*) wavN,wavDT
        allocate(fltr4(wavN),waveTmp(wavN))
        fltr4=0.0
        read(10,*) (fltr4(i),i=1,wavN)
        fltr4=fltr4-fltr4(1)
        
        ! Filter and integrate waveform
        dtr4=real(wavDT)
        if(f1r4>0.)then
          CALL XAPIIR(fltr4, wavN, 'BU', 0.0, 0.0, 4,'BP', f1r4, f2r4, dtr4, 1, wavN)
        else
          CALL XAPIIR(fltr4, wavN, 'BU', 0.0, 0.0, 4,'LP', f1r4, f2r4, dtr4, 1, wavN)
        endif
        
        ! Time integration
        write(*,*) '   -> time integration'
        do mm=2,wavN
          fltr4(mm)=fltr4(mm)+fltr4(mm-1)
        enddo
        fltr4=real(fltr4*wavDT)
        
        ! Time integration (the second just for K-net and KiK-net)
        if(wavDT.le.0.02d0)then
          write(*,*) '   -> time integration'
          do mm=2,wavN
            fltr4(mm)=fltr4(mm)+fltr4(mm-1)
          enddo
          fltr4=fltr4*wavDT
        endif
        
        ! Write waveforms to file for debugging
        if(k==-1)then
          write(jjtxt,'(I0)') jj
          filename='tempWave'//trim(jjtxt)//'.tmp'
          open(230,file=filename)
          do mm=1,wavN
            write(230,*)wavDT*(mm-1),fltr4(mm)
          enddo
          close(230)
        endif
        
        ! Downsample and fill Uvec
        write(*,'(A,F5.1,A)') '    -> downsample (take every',(dt/wavDT),' sample)'
        dspl = int(dt/wavDT)
        waveTmp=0.0
        cou = 0
        do mm=1,wavN
          if(mod(mm,dspl)==1)then
            cou = cou + 1
            waveTmp(cou) = fltr4(mm)
          endif
        enddo
        Uvec((j-1)*nT+1:j*nT) = dble(waveTmp(iT1:iT1+nT-1))
        
        deallocate(fltr4,waveTmp)
        
        enddo
    enddo
    close(10)
    
    ! Write final Unez
    write(*,*) "------------------------------------------------------"
    write(*,*)' Write Uvec data in ',trim(input_Uvec)
    open(12,file=trim(input_Uvec))
    do i=1,nT
      write(12,'(1000E17.7E3)')dt*(iT1-1+i-1),(Uvec((j-1)*nT+i),j=1,NSTAcomp)
    enddo
    close(12)
    
    write(*,*) "DONE"
!---------------------------------------------------------------------
END PROGRAM

