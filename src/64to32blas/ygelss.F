      SUBROUTINE YGELSS( M, N, NRHS, A, LDA,  B, LDB, S, RCOND, RANK,
     W     WORK, LWORK, INFO )
*
* $Id$
*     .. Scalar Arguments ..
      implicit none
#include "y64.fh"
      INTEGER            INFO, LDA, LDB, M, N, NRHS, LWORK, RANK
      double precision   rcond
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK(*), S(*)
c
      INTGR4            INFO4, LDA4, LDB4, M4, N4, NRHS4, LWORK4,
     C     rank4
c
      LDA4=lda
      ldb4=ldb
      lwork4=lwork
      rank4=rank
      m4=m
      n4=n
      nrhs4=nrhs
c
      call DGELSs(m4, N4, NRHS4, A, LDA4, B, LDB4, s, rcond, rank4,
     W     WORK, LWORK4, INFO4 )
c
      info=info4
      RETURN
      END
