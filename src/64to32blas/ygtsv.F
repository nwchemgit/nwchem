      SUBROUTINE YGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
*
* $Id$
*     .. Scalar Arguments ..
      implicit none
#include "y64.fh"
      INTEGER            INFO, LDA, LDB, N, NRHS

      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
c
      INTGR4            INFO4,  LDB4, N4, NRHS4

      ldb4=ldb
      n4=n
      nrhs4=nrhs
c
      call DGTSV( N4, NRHS4, DL, D, DU, B, LDB4, INFO4 )
c
      info=info4
      RETURN
      END
