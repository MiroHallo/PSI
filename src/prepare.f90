!   Read inputs and prepare points on fault-planes for GFs
!   Original code by Frantisek Gallovic
!   Modified by Miroslav Hallo, 2017

PROGRAM PREPARE
!---------------------------------------------------------------------
!  Main program
!  Prepare .dat files with sources on fault planes for GF computation
!---------------------------------------------------------------------
    implicit none

    integer, parameter :: ngmax=10000
    real, parameter :: pi=3.141592654
    real df,aw1,t0
    complex rseis,ui,freq
    real,allocatable,dimension(:):: strike,dip,rake,leng,widt,hhypo
    real,allocatable,dimension(:,:):: hypo,offsets,lrake
    real TM(3,3),ITM(3,3)
    real NEZhypo(3),hypo2(3),xi(3),sour(3)
    real dx1,dx2
    real,allocatable:: x1a(:),x2a(:)
    integer,allocatable:: ng1(:),ng2(:)
    integer NSeg,np,nrake
    integer i,j,k,kk,r
    
    interface
      function Transf(NEZ,smer)
        logical smer
        real Transf(3)
        real NEZ(3)
      end function Transf
    end interface

    common /transform/ TM,NEZhypo,ITM,hypo2
    integer nc,nfreq,nr,ns,ikmax
    real tl,aw,xl,uconv,fref
    namelist  /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,fref
    character(len=6) filename

    !----------------------------------
    ! Read input file
    open(1,file='input/input.dat')

    read(1,*)
    read(1,*) nfreq
    read(1,*)
    read(1,*) np
    read(1,*)
    read(1,*) tl
    read(1,*)
    read(1,*) t0
    read(1,*)
    read(1,*) nr
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*) NSeg
    allocate(ng1(NSeg),ng2(NSeg),strike(NSeg),dip(NSeg),rake(NSeg))
    allocate(hhypo(NSeg),leng(NSeg),widt(NSeg),hypo(3,NSeg))
    allocate(offsets(2,NSeg))
    read(1,*)
    read(1,*) nrake
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*) (strike(kk),dip(kk),rake(kk),kk=1,NSeg)
    read(1,*)
    read(1,*) (leng(kk),widt(kk),kk=1,NSeg)
    read(1,*)
    read(1,*) (ng2(kk),ng1(kk),kk=1,NSeg)
    read(1,*)
    read(1,*) (hypo(2,kk),hypo(1,kk),kk=1,NSeg)
    read(1,*)
    read(1,*) (hhypo(kk),kk=1,NSeg)
    read(1,*)
    read(1,*) (offsets(1,kk),offsets(2,kk),kk=1,NSeg)
    close(1)
    
    !----------------------------------
    ! Process rake-angle scatter
    allocate(lrake(nrake,NSeg))
    if(nrake==1)then
      write(*,*) '  -> Fixed rake'
      do kk=1,NSeg
        lrake(1,kk)=rake(kk)
      enddo
    elseif(nrake==2)then
      write(*,*) '  -> Allowed rake scatter (two orthogonal rakes)'
      do kk=1,NSeg
        lrake(1,kk)=rake(kk)-45.
        if(lrake(1,kk)<-180.)then
          lrake(1,kk)=lrake(1,kk)+360.
        endif
        lrake(2,kk)=rake(kk)+45.
        if(lrake(2,kk)>180.)then
          lrake(2,kk)=lrake(2,kk)-360.
        endif
      enddo
    else
      write(*,*) ' Error! Not supported number of rakes!!!'
    endif
    
    !----------------------------------
    ! Prepare discrete points on faults
    open(1,file='sources.dat')
    open(2,file='fault.dat')
    
    k=0
    do r=1,nrake
      do kk=1,NSeg
        if(ng2(kk)>ngmax.or.ng1(kk)>ngmax)stop 'Check dimensions!'

        hypo(3,kk)=0.

        allocate(x1a(ng1(kk)),x2a(ng2(kk)))

        NEZhypo(1)=offsets(1,kk)
        NEZhypo(2)=offsets(2,kk)
        NEZhypo(3)=-hhypo(kk)
        hypo2=hypo(:,kk)
      
        TM(1,1)=sin(strike(kk)/180.0*pi)*cos(dip(kk)/180.0*pi)
        TM(1,2)=-cos(strike(kk)/180.0*pi)*cos(dip(kk)/180.0*pi)
        TM(1,3)=sin(dip(kk)/180.0*pi)
        TM(2,1)=cos(strike(kk)/180.0*pi)
        TM(2,2)=sin(strike(kk)/180.0*pi)
        TM(2,3)=0.
        TM(3,1)=-sin(strike(kk)/180.0*pi)*sin(dip(kk)/180.0*pi)
        TM(3,2)=cos(strike(kk)/180.0*pi)*sin(dip(kk)/180.0*pi)
        TM(3,3)=cos(dip(kk)/180.0*pi)
        ITM=transpose(TM)

        dx1=widt(kk)/float(ng1(kk))
        dx2=leng(kk)/float(ng2(kk))

        do i=1,ng1(kk)
          x1a(i)=(float(i)-.5)*dx1
        enddo
        do i=1,ng2(kk)
          x2a(i)=(float(i)-.5)*dx2
        enddo
      
        do i=1,ng1(kk)
          do j=1,ng2(kk)
            xi(1)=x1a(i)
            xi(2)=x2a(j)
            xi(3)=0.
            k=k+1
            if(k<10)then
              write(filename,'(A5,I1)')'00000',k
            elseif(k<100)then
              write(filename,'(A4,I2)')'0000',k
            elseif(k<1000)then
              write(filename,'(A3,I3)')'000',k
            elseif(k<10000)then
               write(filename,'(A2,I4)')'00',k
            elseif(k<100000)then
              write(filename,'(A1,I5)')'0',k
            else
              write(filename,'(I6)')k
            endif
            sour=Transf(xi,.FALSE.)
            sour(3)=-sour(3)
            sour=sour/1000.
            write(1,'(A6,6E13.5)') filename,sour,strike(kk),dip(kk),lrake(r,kk)
          enddo
        enddo
      
        deallocate(x1a,x2a)
        
        if(r>1)then
          cycle
        endif
        
        xi=0.
        sour=Transf(xi,.FALSE.)/1000.
        sour(3)=-sour(3)
        write(2,'(3F11.4)') sour
        xi(1)=widt(kk)
        xi(2)=0.
        xi(3)=0.
        sour=Transf(xi,.FALSE.)/1000.
        sour(3)=-sour(3)
        write(2,'(3F11.4)') sour
        xi(1)=widt(kk)
        xi(2)=leng(kk)
        xi(3)=0.
        sour=Transf(xi,.FALSE.)/1000.
        sour(3)=-sour(3)
        write(2,'(3F11.4)') sour
        xi(1)=0.
        xi(2)=leng(kk)
        xi(3)=0.
        sour=Transf(xi,.FALSE.)/1000.
        sour(3)=-sour(3)
        write(2,'(3F11.4)') sour
        xi=0.
        sour=Transf(xi,.FALSE.)/1000.
        sour(3)=-sour(3)
        write(2,'(3F11.4)') sour
        write(2,*)
      
      enddo
    enddo
    
    close(1) 
    close(2)
    
    !----------------------------------
    ! Prepare help file for dwn

    open(2,file='GRDAT.HED')
    aw=1.;ns=1;xl=3000000.;ikmax=200000;uconv=1.E-12;fref=1. !Axitra values that does not have to be generally changed
    nc=0;  !set up formal value that are actualy not used by Axitra (they are readed from elsewhere)
    write(2,input)
    close(2)

    ui=cmplx(0.,1.)
    df=1./tl
    aw1=-aw/(2.*tl)
    open(1,file='dirac.dat')
    do i=1,nfreq
      freq=cmplx(df*(i-1),aw1)
      rseis=exp(-ui*2.*pi*t0*freq)
      write(1,*) real(rseis),imag(rseis)
    enddo
    do i=nfreq+1,np
      rseis=0.
      write(1,*) real(rseis),imag(rseis)
    enddo
    close(1)

    open(3,file='sortseg.dat')
    write(3,*) "No. of computed frequencies"
    write(3,'(I6)') nfreq
    write(3,*) "Number of segments"
    write(3,'(I4)') NSeg
    write(3,*) "Number of receivers"
    write(3,'(I4)') nr
    write(3,*) "Number of rakes"
    write(3,'(I4)') nrake
    write(3,*) "Computed points along strike and dip"
    write(3,'(2I5)') (ng2(kk),ng1(kk),kk=1,NSeg)
    write(3,*) "Strike  Dip   Rake"
    write(3,'(3F9.2)') (strike(kk),dip(kk),rake(kk),kk=1,NSeg)
    write(3,*) "Length and width"
    write(3,'(2F11.2)') (leng(kk),widt(kk),kk=1,NSeg)
    write(3,*) "Position of reference point on the fault"
    write(3,'(2F11.2)') (hypo(2,kk),hypo(1,kk),kk=1,NSeg)
    write(3,*) "Depth of fault reference point"
    write(3,'(F11.2)') (hhypo(kk),kk=1,NSeg)
    close(3)

    deallocate(strike,dip,rake,leng,widt,hhypo,hypo,offsets,ng1,ng2)
    deallocate(lrake)
!---------------------------------------------------------------------
END PROGRAM




FUNCTION Transf(NEZ,smer)
!---------------------------------------------------------------------
    implicit none
    logical smer
    real Transf(3)
    real NEZhypo(3),hypo2(3)
    real TM(3,3),ITM(3,3)
    real NEZ(3)
    common /transform/ TM,NEZhypo,ITM,hypo2

    if (smer) then
      Transf=matmul(TM,(NEZ-NEZhypo))+hypo2
    else
      Transf=matmul(ITM,(NEZ-hypo2))+NEZhypo
    endif
!---------------------------------------------------------------------
END FUNCTION
