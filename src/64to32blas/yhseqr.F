      SUBROUTINE YHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
     $                   LDZ, WORK, LWORK, INFO )
c $Id$
*
      implicit none
#include "y64.fh"
      CHARACTER          COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   Z( LDZ, * )
c
      INTGR4            IHI4, ILO4, INFO4, LDH4, LDZ4, LWORK4, N4
c
      IHI4=ihi
      ILO4=ilo
      LDH4=ldh
      LDZ4=ldz
      LWORK4=lwork
      N4=n
c
      call DHSEQR( JOB, COMPZ, N4, ILO4, IHI4, H, LDH4, WR, WI, Z,
     $                   LDZ4, WORK, LWORK4, INFO4 )
      info=info4
      return
      END
