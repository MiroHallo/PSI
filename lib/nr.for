      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
      
      
      
!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*        N      size of table RA                   *
!*          RA      table to be sorted                 *
!* OUTPUT:                                           *
!*        RA    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
      SUBROUTINE HPSORT(N,RA)
      Double precision RA(N)
      Double precision RRA
      L=N/2+1
      IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR 
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
   10 continue
      if(L > 1)then
      L=L-1
      RRA=RA(L)
      else
      RRA=RA(IR)
      RA(IR)=RA(1)
      IR=IR-1
      if(IR.eq.1)then
      RA(1)=RRA
      return
      end if
      end if
      I=L
      J=L+L
   20 if(J.le.IR)then
      if(J < IR)then
      if(RA(J) < RA(J+1))  J=J+1
      end if
      if(RRA < RA(J))then
      RA(I)=RA(J)
      I=J; J=J+J
      else
      J=IR+1
      end if
      goto 20
      end if
      RA(I)=RRA
      goto 10
      END
      
      
      
!-------------------------------------------------------------------
!                                               
!       Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ------------------------------------------------------------------

      FUNCTION GASDEV(idum)

!     ..Arguments..
      integer          idum
      real GASDEV

!     ..Local..
      real v1,v2,r,fac
      real ran3

      if (idum.lt.0) iset=0
 10   v1=2*ran3(idum)-1
      v2=2*ran3(idum)-1
      r=v1**2+v2**2
      if(r.ge.1.or.r.eq.0) GOTO 10
      fac=sqrt(-2*log(r)/r)
      GASDEV=v2*fac

      RETURN
      END
      
      
      
!-------------------------------------------------------------------
!                                               
!   Matrix inverse using Cholesky decomposition
!
!   input    n  size of matrix
!   input     A  Symmetric positive def. matrix
!   output  aa  inverse of A
!   uses        choldc1(int,MAT,VEC)
!-------------------------------------------------------------------

      SUBROUTINE cholsl(n,A,aa)
      integer n
      real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
      integer i,j,k
      
      call choldcsl(n,A,aa)
      do i = 0, n-1
        do j = i + 1, n-1
          aa(i,j) = 0.d0
        enddo
      enddo
      do i = 0, n-1
        aa(i,i) = aa(i,i) * aa(i,i)
        do k = i + 1, n-1
          aa(i,i) = aa(i,i) + aa(k,i) * aa(k,i)
        enddo
        do j = i + 1, n-1
          do k = j, n-1
            aa(i,j) = aa(i,j) + aa(k,i) * aa(k,j)
          enddo
        enddo
      enddo
      do i = 0,  n-1
        do j = 0, i-1
          aa(i,j) = aa(j,i)
        enddo
      enddo
      return
      END
      
      
      
! -----------------------------------------------------
!     Inverse of Cholesky decomposition.
!
!     input    n  size of matrix
!     input    A  Symmetric positive def. matrix
!     output  aa  inverse of lower decomposed matrix
!     uses        choldc1(int,MAT,VEC)         
! -----------------------------------------------------
      
      SUBROUTINE choldcsl(n,A,aa)
      integer n
      real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
      integer i,j,k, ialloc
      real*8 sum
      real*8, pointer :: p(:)
      allocate(p(0:n-1),stat=ialloc)
      
      aa = A
      call choldc1(n, aa, p)
      do i = 0, n-1
        aa(i,i) = 1.d0 / p(i)
        do j = i + 1, n-1
          sum = 0.d0
          do k = i, j-1
            sum = sum - aa(j,k) * aa(k,i)
          enddo
          aa(j,i) = sum / p(j)
        enddo
      enddo
      deallocate(p)
      return
      END
      
      
      
! -------------------------------------------------
! main method for Cholesky decomposition.
!
! input         n  size of matrix
! input/output  a  matrix
! output        p  vector of resulting diag of a
! -------------------------------------------------

      SUBROUTINE choldc1(n,a,p)
      integer n
      real*8 a(0:n-1,0:n-1), p(0:n-1)
      integer i,j,k
      real*8 sum
      
      do i = 0, n-1
        do j = i, n-1
          sum = a(i,j)
          do k = i - 1, 0, -1
            sum = sum - a(i,k) * a(j,k)
          enddo
          if (i.eq.j) then
            if (sum <= 0.d0) then
              write(*,*)sum
              print *,' Error! the matrix is not positive definite!'
            endif
            p(i) = dsqrt(sum)
          else
            a(j,i) = sum / p(i)
          endif
        enddo
      enddo
      return
      END