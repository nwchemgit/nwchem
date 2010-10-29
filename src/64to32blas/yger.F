      SUBROUTINE YGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
* $Id$
*
      implicit none
#include "y64.fh"
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
c
      INTGR4          INCX4, INCY4, LDA4, M4, N4
c
      INCX4=incx 
      INCY4=incy 
      LDA4=lda 
      M4=m 
      N4=n
c
      call DGER(M4, N4, ALPHA, X, INCX4, Y, INCY4, A, LDA4 )
c
      RETURN
      END
