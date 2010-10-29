      SUBROUTINE YSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), B(LDB,*), W( * ), WORK( * )
c
      INTGR4            INFO4, LDA4, LDB4, LWORK4, N4,ITYPE4
c
      lda4=lda
      ldb4=ldb
      lwork4=lwork
      itype4=itype
      n4=n
      call DSYGV( itype4, JOBZ, UPLO, N4, A, LDA4, B, ldb4, W, WORK, 
     $     LWORK4, INFO4)
      info=info4
      RETURN
      END
