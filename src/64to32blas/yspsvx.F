      SUBROUTINE YSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X,
     $                   LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          FACT, UPLO
      DOUBLE PRECISION   RCOND
      INTEGER            INFO, LDB, LDX, N, NRHS
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   AFP( * ), AP( * ), B( LDB, * ), BERR( * ),
     $                   FERR( * ), WORK( * ), X( LDX, * )
c
      INTGR4            INFO4, LDB4, LDX4, N4, NRHS4
C Warning the values of IPIV should not be used, since they are not
C changed before outputting
cedo      INTEGER            IPIV4(maxn), IWORK4(maxn)
c
      ldb4=ldb
      ldx4=ldx
      n4=n
      nrhs4=nrhs
c
      call DSPSVX( FACT, UPLO, N4, NRHS4, AP, AFP, IPIV, B, LDB4, X,
     $                   LDX4, RCOND, FERR, BERR, WORK, IWORK, INFO4)
      info=info4
c
      RETURN
      END
