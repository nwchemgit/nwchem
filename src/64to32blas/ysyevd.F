      SUBROUTINE YSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
     *     LIWORK, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, LIWORK,N
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
      INTEGER            IWORK(*)
c
      INTGR4            INFO4, LDA4, LWORK4, N4, LIWORK4
c
      lda4=lda
      lwork4=lwork
      liwork4=liwork
      n4=n
      call DSYEVD( JOBZ, UPLO, N4, A, LDA4, W, WORK, LWORK4, IWORK,
     *     LIWORK4, INFO4)
      info=info4
      RETURN
      END
