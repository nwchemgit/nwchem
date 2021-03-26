      SUBROUTINE YPOTRF( UPLO, N, A, LDA, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c
      INTGR4 n4,lda4,info4
c
      n4=n
      lda4=lda
c
      call DPOTRF( UPLO, N4, A, LDA4, INFO4 )
      info=info4
      return
      end
