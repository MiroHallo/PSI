!---------------------------------------------------------------------
!
!  Creates synthetic slip-rate time source functions (tfin)
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

SUBROUTINE synth_tfin2()
!---------------------------------------------------------------------
!  Creates synthetic slip-rate time source functions by splines
!  Result in tfin(:) in CoreBox and save to files
!  if syntdata==2
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use ExeBox
    implicit none
    character(len=80):: fileprefix
    integer im,kk,error
    
    !----------------------------------
    !  Allocate
    im=1
    allocate(splineN0(NSeg))
    do kk=1,NSeg
      splineN0(kk)=2*NL(kk)+2*NW(kk)
    enddo
    allocate(m_veloc(maxval(NPL),maxval(NPW),NSeg,im))
    allocate(m_risetime(maxval(NPL),maxval(NPW),NSeg,im))
    allocate(m_peaktime(maxval(NPL),maxval(NPW),NSeg,im))
    allocate(m_rakes(maxval(NPL),maxval(NPW),NSeg,im))
    allocate(m_hypo(3,NSeg,im))
    allocate(m_slip(maxval(splineN0)+MaxNSpline,3,NSeg,im))
    allocate(m_Nspline(NSeg,im))
    allocate(slip(maxval(NL),maxval(NW),NSeg))
    allocate(risetime(maxval(NL),maxval(NW),NSeg))
    allocate(ruptime(maxval(NL),maxval(NW),NSeg))
    allocate(peaktime(maxval(NL),maxval(NW),NSeg))
    allocate(rakes(maxval(NL),maxval(NW),NSeg))
    allocate(timefunc(Ssvd,maxval(NL),maxval(NW),NSeg))
    
    !----------------------------------
    ! Constant model parameters
    m_veloc = 3000.d0
    m_risetime = 2.6d0
    m_peaktime = 0.8d0
    do kk=1,NSeg
      m_rakes(:,:,kk,im) = rake(kk) -20
    enddo
    
    !----------------------------------
    ! Hypocentrum (rupture start point)
    do kk=1,NSeg
      m_hypo(1,kk,im) = epicL(kk)
      m_hypo(2,kk,im) = epicW(kk)
      m_hypo(3,kk,im) = segDelay(kk)
    enddo
    
    !----------------------------------
    ! Slip by spline model parameters
    m_slip=0.d0
    do kk=1,NSeg
      ! Spline zero control points on edges
      m_slip(1:NL(kk),1,kk,im)=cellL(1:NL(kk),kk)
      m_slip(NL(kk)+1:2*NL(kk),1,kk,im)=cellL(1:NL(kk),kk)
      m_slip(NL(kk)+1:2*NL(kk),2,kk,im)=widt(kk)
      m_slip(2*NL(kk)+1:2*NL(kk)+NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
      m_slip(2*NL(kk)+NW(kk)+1:2*NL(kk)+2*NW(kk),1,kk,im)=leng(kk)
      m_slip(2*NL(kk)+NW(kk)+1:2*NL(kk)+2*NW(kk),2,kk,im)=cellW(1:NW(kk),kk)
      if(kk.eq.1)then
        ! Three spline efficient control points
        m_Nspline(kk,im) = 3
        m_slip(splineN0(kk)+1,1,kk,im) = leng(kk)/4.d0
        m_slip(splineN0(kk)+1,2,kk,im) = widt(kk)/2.d0
        m_slip(splineN0(kk)+1,3,kk,im) = 1.d0
        m_slip(splineN0(kk)+2,1,kk,im) = leng(kk)/2.d0
        m_slip(splineN0(kk)+2,2,kk,im) = widt(kk)/2.d0
        m_slip(splineN0(kk)+2,3,kk,im) = -1.d0
        m_slip(splineN0(kk)+3,1,kk,im) = leng(kk)/1.333d0
        m_slip(splineN0(kk)+3,2,kk,im) = widt(kk)/2.d0
        m_slip(splineN0(kk)+3,3,kk,im) = 1.d0
      else
        ! One spline efficient control point
        m_Nspline(kk,im) = 1
        m_slip(splineN0(kk)+1,1,kk,im) = leng(kk)/2.d0
        m_slip(splineN0(kk)+1,2,kk,im) = widt(kk)/2.d0
        m_slip(splineN0(kk)+1,3,kk,im) = 0.5d0
      endif
    enddo
    
    !----------------------------------
    ! Interpret model parameters
    call InterpMod(im,Mfix,error)
    
    !----------------------------------
    !  If error
    if(error.ne.0)then
      write(*,*)' Error! during interpretation model parameters; GenErr:',error
      stop
    endif
    
    !----------------------------------
    !  Save slip-rate model
    write(*,'(A,E10.2,A)')' Save slip-rate model (M0:',moment,')'
    fileprefix = 'inv/synth'
    call SaveSlipMod(fileprefix)
    
    !----------------------------------
    !  Deallocate
    deallocate(m_veloc,m_risetime,m_peaktime,m_rakes,m_hypo,m_slip,m_Nspline,splineN0)
    deallocate(timefunc,slip,risetime,ruptime,peaktime,rakes)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE synth_tfin1()
!---------------------------------------------------------------------
!  Creates synthetic slip-rate time source functions
!  Result in tfin(:) in CoreBox and save to files
!  if syntdata==1
!---------------------------------------------------------------------
    use InputBox
    use CoreBox
    use LibBox
    use ExeBox
    implicit none
    integer,allocatable,dimension(:):: syntdatai,syntdataj
    real(8),allocatable,dimension(:):: syntdelay
    character(len=80):: fileprefix
    real(8) posL,posW,hypL,hypW,dum,mome,rangd
    integer i,j,k,kk,r,coun,Stype
    
    !----------------------------------
    !  Allocate
    allocate(syntdatai(NSeg),syntdataj(NSeg),syntdelay(NSeg))
    allocate(slip(maxval(NL),maxval(NW),NSeg))
    allocate(risetime(maxval(NL),maxval(NW),NSeg))
    allocate(ruptime(maxval(NL),maxval(NW),NSeg))
    allocate(peaktime(maxval(NL),maxval(NW),NSeg))
    allocate(rakes(maxval(NL),maxval(NW),NSeg))
    allocate(timefunc(Ssvd,maxval(NL),maxval(NW),NSeg))
    
    !----------------------------------
    !  SET INPUT
    Stype = 3 ! 1=delta, 2=single asper, 3= double asper.
    syntdatai = (/ 0, 0 /) ! Point along strike
    syntdataj = (/ 0, 0 /) ! Point along dip (from down to up)
    syntdelay = (/ 0.d0, 1.5d0 /) ! Segment time delay [sec] for synthetics
    
    risetime = 2.6d0
    peaktime = 0.8d0
    
    !----------------------------------
    !  Delta function (slip and ruptime)
    if(Stype==1)then
      write(*,*)' Delta function slip-rate for resolution anal.'
      slip=0.d0
      if(syntdatai(1)==0.or.syntdataj(1)==0)then
        write(*,*)' Error! Undefined position of patch with slip!'
        stop
      elseif(syntdatai(1)>NL(1).or.syntdataj(1)>NW(1))then
        write(*,*)' Error! Hypocentre position out of range!'
        stop
      endif
      slip(syntdatai(1),syntdataj(1),1) = 1.d0
      ruptime(syntdatai(1),syntdataj(1),1) = 2.d0
    
    !----------------------------------
    !  Full slip-rate model (slip and ruptime)
    else
      write(*,*)' Synthetic full slip-rate model'
      slip=1.d0
      do kk=1,NSeg
      
        ! Set up the hypocenter position (from reference point or manual)
        if(syntdatai(kk)==0.or.syntdataj(kk)==0)then
          hypL=epicL(kk)
          hypW=epicW(kk)
        elseif(syntdatai(kk)>NL(kk).or.syntdataj(kk)>NW(kk))then
          write(*,*)'  --> Error! Hypocentre position out of range!'
          stop
        else
          hypL=(dble(syntdatai(kk))-0.5d0)*leng(kk)/dble(NL(kk))
          hypW=(dble(syntdataj(kk))-0.5d0)*widt(kk)/dble(NW(kk))
        endif
        
        ! Prepare ruptime
        do j=1,NW(kk)
          posW=(dble(j)-0.5d0)*widt(kk)/dble(NW(kk))
          do i=1,NL(kk)
            posL=(dble(i)-0.5d0)*leng(kk)/dble(NL(kk))
            ruptime(i,j,kk)=syntdelay(kk)+sqrt((posL-hypL)**2+(posW-hypW)**2)/3000
            if(ruptime(i,j,kk)+risetime(i,j,kk)+1.575*peaktime(i,j,kk)>dble(Ssvd-1)*dt)then
              write(*,*)' Warning! Source time function longer than time window!'
            endif
          enddo
        enddo
        
        ! Prepare slip
        do j=1,NW(kk)
          posW=(dble(j)-0.5d0)*widt(kk)/dble(NW(kk))
          do i=1,NL(kk)
            posL=(dble(i)-0.5d0)*leng(kk)/dble(NL(kk))
            if(Stype==2)then ! single elliptical crack
              slip(i,j,kk)=exp(-(posL-leng(kk)/2)**2/(leng(kk)/3)**2) * exp(-(posW-widt(kk)/2)**2/(widt(kk)/3)**2)
            elseif(Stype==3)then ! two asperities in the first segment, single elliptical cracks in next segments
              if(kk==1)then
                slip(i,j,kk)=exp(-(posL-leng(kk)/4)**2/(leng(kk)/7)**2) * exp(-(posW-widt(kk)/2)**2/(widt(kk)/4)**2) + &
                           exp(-(posL-leng(kk)*3/4)**2/(leng(kk)/7)**2) * exp(-(posW-widt(kk)/2)**2/(widt(kk)/4)**2)
              else
                slip(i,j,kk)=exp(-(posL-leng(kk)/2)**2/(leng(kk)/3)**2) * exp(-(posW-widt(kk)/2)**2/(widt(kk)/3)**2)
              endif
            else
              write(*,*)' Error! Undefined type of slip-rate model!!!'
              stop
            endif
          enddo
        enddo
        
      enddo
    endif
    
    !----------------------------------
    !  Prepare rakes
    do kk=1,NSeg
      do j=1,NW(kk)
        posW=(dble(j)-0.5d0)*widt(kk)/dble(NW(kk))
        do i=1,NL(kk)
          posL=(dble(i)-0.5d0)*leng(kk)/dble(NL(kk))
          rakes(i,j,kk)=rake(kk)
        enddo
      enddo
    enddo
    
    !----------------------------------
    !  Generates slip-rate time function
    timefunc=0.d0
    do kk=1,NSeg
      do j=1,NW(kk)
        posW=(dble(j)-.5d0)*dW(kk)
        do i=1,NL(kk)
          posL=(dble(i)-.5d0)*dL(kk)
          ! Regularized Yoffe slip-rate time function
          timefunc(1:Ssvd,i,j,kk)=RYoffe(Ssvd,dt,ruptime(i,j,kk),risetime(i,j,kk),peaktime(i,j,kk),slip(i,j,kk))
          ! In case of GPS inversion only
          if(Ssvd==1)timefunc(1,i,j,kk)=slip(i,j,kk)    
          ! Set ruptime to 0 for pathes without STF
          dum=sum(real(timefunc(:,i,j,kk)))*dt
          if(dum==0.d0)ruptime(i,j,kk)=0.d0
        enddo
      enddo
    enddo
    
    ! Norm. total seismic moment
    mome=0.d0
    do kk=1,NSeg
      do k=1,Ssvd
        mome=mome+sum(timefunc(k,1:NL(kk),1:NW(kk),kk)*mu(1:NL(kk),1:NW(kk),kk))*elem(kk)*dt
      enddo
    enddo
    timefunc=timefunc/mome*Mfix
    
    !----------------------------------
    !  Regroup timefunc into tfin and apply rake coef.
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
    
    !----------------------------------
    !  Save slip-rate model
    write(*,*)' Save slip-rate model'
    fileprefix = 'inv/synth'
    call SaveSlipMod(fileprefix)
    
    !----------------------------------
    !  Deallocate
    deallocate(timefunc,slip,risetime,ruptime,peaktime,rakes)
    deallocate(syntdatai,syntdataj,syntdelay)
    
!---------------------------------------------------------------------
END SUBROUTINE


