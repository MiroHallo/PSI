!   Read results from dwn and combine all point into one file
!   Original code by Frantisek Gallovic
!   Modified by Miroslav Hallo, 2017

PROGRAM  RESORT
!---------------------------------------------------------------------
!  Main program
!  Prepare NEZsor.dat file with all elem. seismograms from all points
!---------------------------------------------------------------------
	implicit none

	integer i,k,kk,r
	integer np,nr,rec,nfmax,NSeg,gntot,grtot,nrake
	integer,allocatable:: gn1(:),gn2(:)
	real*4 dat(7)
	CHARACTER*6 filename
	
	!----------------------------------
	! Read input/output files
	open(3,form='unformatted',file='NEZsor.dat')
	open(5,file='sortseg.dat',action='read',status='old')

	read(5,*)
	read(5,*) nfmax
	read(5,*)
	read(5,*) NSeg
	allocate(gn1(NSeg),gn2(NSeg))
	read(5,*)
	read(5,*) nr
	read(5,*)
	read(5,*) nrake
	read(5,*)
	read(5,*) (gn2(kk),gn1(kk),kk=1,NSeg)
	close(5)
	
	!----------------------------------
	! Index all files with results
	gntot=sum(gn1(:)*gn2(:))*nrake
	do i=1,gntot
	  if(i<10)then
		write(filename,'(A5,I1)')'00000',i
	  elseif(i<100)then
		write(filename,'(A4,I2)')'0000',i
	  elseif(i<1000)then
		write(filename,'(A3,I3)')'000',i
	  elseif(i<10000)then
		write(filename,'(A2,I4)')'00',i
	  elseif(i<100000)then
		write(filename,'(A1,I5)')'0',i
	  else
		write(filename,'(I6)')i
	  endif
	  open(100+i,FILE='dat/'//filename//'.nez',action='read',status='old')
	enddo
	
	!----------------------------------
	! Write into one file
	write(*,*) "No.frequencies, No.receivers, No.segments, No.rakes:"
	write(*,*) nfmax,nr,NSeg,nrake
	grtot=sum(gn1(:)*gn2(:))
	do r=1,nrake
	  write(*,*) "Rake No. ",r
	  do kk=1,NSeg
	    write(*,*) "Segment ",kk
	    do rec=1,nr
		  write(*,*) rec
		  write(3) rec
		  do i=((r-1)*grtot)+sum(gn1(1:kk-1)*gn2(1:kk-1))+1,((r-1)*grtot)+sum(gn1(1:kk)*gn2(1:kk))
		    read(100+i,*)
		    write(3) i
		    do k=1,nfmax
			  read(100+i,*) dat(1:7)
			  write(3) dat(2:7)
		    enddo
		  enddo
	    enddo
	  enddo
	enddo

	close(3)

END PROGRAM




