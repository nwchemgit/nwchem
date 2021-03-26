      subroutine ysymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
*
* $Id$
*
      implicit none
#include "y64.fh"
      double precision alpha,beta
      integer incx,incy,lda,n
      character uplo
      double precision a(lda,*),x(*),y(*)
c
      INTGR4 n4,lda4,incx4,incy4
c
      n4=n
      lda4=lda
      incx4=incx
      incy4=incy
      call dsymv(uplo,n4,alpha,a,lda4,x,incx4,beta,y,incy4)
      return
      end
