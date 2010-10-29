*
* $Id$
*
      subroutine types1(a,b,alpha,beta,q,e,NIN,NDIM)
C
C-----------------------------------------------------------------------
C
C     This routine reads Fann's matrices, i.e., alpha's/beta's,
C     and computes the representation a's/b's, and q's/e's. 
C
C-----------------------------------------------------------------------
C
      integer          NDIM,NIN
      double precision a(*),alpha(*),b(*),beta(*),e(*),q(*)
C
      double precision zero
      parameter        (zero=0.0d0)
      integer          i,j
C
      intrinsic        sqrt
C
C.... read the alpha's and beta's ......................................
C
cgif      call parser (record,string)
cgif      read (string,fmt='(a)') matrix
cgif      call parser (record,string)
cgif      read (string,fmt='(e24.0)') shift
C
cgif      open (9,err=1,file=matrix,form='formatted')
cgif      read (9,err=1,fmt=*) NIN 
cgif      if ( NIN.GT.NDIM ) stop 'TYPE1: not enough space in Z'
cgif      read (9,err=1,fmt=*) (j,alpha(i),beta(i),i=1,NIN)
C
cgif      do j = 1,NIN-1
cgif         alpha(j) = alpha(j) + shift
cgif         beta(j) = beta(j+1)
cgif      end do
cgif      alpha(NIN) = alpha(NIN) + shift
C
C.... compute the q's and e's ..........................................
C
      q(1) = alpha(1)
      do j = 1,NIN-1
         e(j) = (beta(j)/q(j)) * beta(j)
         q(j+1) = alpha(j+1) - e(j) 
         if ( q(j+1).lt.zero ) 
     &      STOP 'Matrix is not positive definite, change shift'
      end do 
C
C.... compute the a's and b's ..........................................
C
      a(1) = sqrt(q(1))
      do i=1,NIN-1
         b(i) = sqrt(e(i))
         a(i+1) = sqrt(q(i+1))
      end do
C
      return
    1 stop '* I/O error, check file for Fann''s matrix'
      end
