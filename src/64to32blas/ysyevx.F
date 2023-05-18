      SUBROUTINE YSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     *                     ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     *                    IFAIL, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
      
      INTEGER            IFAIL(*), IWORK(*)
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
c
      INTGR4            IL4, INFO4, IU4, LDA4, LDZ4, LWORK4
      INTGR4            M4, N4
c
      il4=il
      iu4=iu
      lda4=lda
      ldz4=ldz
      lwork4=lwork
      m4=m
      n4=n
      CALL DSYEVX( JOBZ, RANGE, UPLO, N4, A, LDA4, VL, VU, IL4, IU4,
     *                     ABSTOL, M4, W, Z, LDZ4, WORK, LWORK4, IWORK,
     *                    IFAIL, INFO4 )
      info=info4
      RETURN
      END
