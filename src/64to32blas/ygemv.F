      SUBROUTINE YGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
* $Id: ygemv.F,v 1.1 2005-11-02 17:14:00 edo Exp $
*
      implicit none
#include "y64.fh"
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
c
      INTGR4            INCX4, INCY4, LDA4, M4, N4
c
      INCX4=incx
      INCY4=incy
      LDA4=lda
      M4=m
      N4=n
c
      call DGEMV ( TRANS, M4, N4, ALPHA, A, LDA4, X, INCX4,
     $                   BETA, Y, INCY4 )

      RETURN
      END
