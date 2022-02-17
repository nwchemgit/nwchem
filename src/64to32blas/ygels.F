      SUBROUTINE YGELS( TRANS, M, N, NRHS, A, LDA,  B, LDB,
     W     WORK, LWORK, INFO )
*
* $Id$
*     .. Scalar Arguments ..
      implicit none
#include "y64.fh"
      character         trans
      INTEGER            INFO, LDA, LDB, M, N, NRHS, LWORK
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK(*)
c
      INTGR4            INFO4, LDA4, LDB4, M4, N4, NRHS4, LWORK4
      integer i
c
      LDA4=lda
      ldb4=ldb
      lwork4=lwork
      m4=m
      n4=n
      nrhs4=nrhs
c
      call DGELS(trans, m4, N4, NRHS4, A, LDA4, B, LDB4,
     W     WORK, LWORK4, INFO4 )
c
      info=info4
      RETURN
      END
