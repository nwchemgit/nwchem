*
* $Id$
*

      SUBROUTINE YPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      implicit none
#include "y64.fh"
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c
      INTGR4            INFO4, LDA4, LDB4, N4, NRHS4
c
      lda4=lda
      ldb4=ldb
      n4=n
      nrhs4=nrhs
c
      call DPOSV( UPLO, N4, NRHS4, A, LDA4, B, LDB4, INFO4)
      info=info4
      RETURN
      END
