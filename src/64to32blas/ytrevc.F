      SUBROUTINE YTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, MM, M, WORK, INFO )
c $Id$
*
      implicit none
#include "y64.fh"
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
C Warning size of logical=size of integer:  Will not work with HOWMNY="S"
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
c
      INTGR4            INFO4, LDT4, LDVL4, LDVR4, M4, MM4, N4
c
      LDT4=ldt
      LDVL4=ldvl
      LDVR4=ldvr
      M4=m
      MM4=mm
      N4=n
c      
      call DTREVC( SIDE, HOWMNY, SELECT, N4, T, LDT4, VL, LDVL4, VR,
     $                   LDVR4, MM4, M4, WORK, INFO4 )
      info=info4
      return
      end
