      SUBROUTINE YSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
c
      INTGR4            INFO4, LDA4, LWORK4, N4
c
      lda4=lda
      lwork4=lwork
      n4=n
      call DSYEV( JOBZ, UPLO, N4, A, LDA4, W, WORK, LWORK4, INFO4)
      info=info4
      RETURN
      END
