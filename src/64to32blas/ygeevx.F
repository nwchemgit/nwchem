      SUBROUTINE YGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
     $                   VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
*
* $Id$
*
      implicit none
#include "y64.fh"
*     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),
     $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
c
      INTGR4            N4, LDA4, LDVL4, LDVR4, LWORK4, INFO4,
     $     ILO4,IHI4
c
      N4=N
      LDA4=LDA
      LDVL4=LDVL
      LDVR4=LDVR
      LWORK4=LWORK
c
      CALL DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N4, A, LDA4, WR, WI,
     $                   VL, LDVL4, VR, LDVR4, ILO4, IHI4, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK4, IWORK, INFO4 )
      INFO=INFO4
      ILO=ILO4
      IHI=IHI4

      RETURN
      END
