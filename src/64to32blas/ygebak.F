      SUBROUTINE YGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
     $                   INFO )
c $Id$
*
      implicit none
#include "y64.fh"
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
c
      INTGR4            IHI4, ILO4, INFO4, LDV4, M4, N4
c
      IHI4=ihi
      ILO4=ilo 
      LDV4=ldv
      M4=m
      N4=n
c
      call DGEBAK( JOB, SIDE, N4, ILO4, IHI4, SCALE, M4, V, LDV4,
     $                   INFO4 )
      info=info4
      return
      END
