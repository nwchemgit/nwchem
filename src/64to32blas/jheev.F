      SUBROUTINE jHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
c
* $Id: xheev.F 19697 2010-10-29 16:57:34Z d3y133 $
c
      implicit none
#include "y64.fh"
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
c
      INTGR4            INFO4, LDA4, LWORK4, N4
c
      lda4=lda
      lwork4=lwork
      n4=n
c
      call ZHEEV( JOBZ, UPLO, N4, A, LDA4, W, WORK, LWORK4, RWORK,
     $                  INFO4 )
      info=info4
      RETURN
      END
