      SUBROUTINE YGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
*
* $Id$
*
      implicit none
#include "y64.fh"
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
c
      INTGR4            N4, LDA4, LDVL4, LDVR4, LWORK4, INFO4
c
      N4=N
      LDA4=LDA
      LDVL4=LDVL
      LDVR4=LDVR
      LWORK4=LWORK
c
      CALL DGEEV( JOBVL, JOBVR, N4, A, LDA4, WR, WI, VL, LDVL4, VR,
     $                  LDVR4, WORK, LWORK4, INFO4 )
      INFO=INFO4

      RETURN
      END
