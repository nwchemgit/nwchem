      Subroutine DEBTest(RTDB)
C$Id: debtest.F,v 1.4 1996-01-17 21:46:38 d3e129 Exp $
      Implicit NONE
      Integer RTDB
      Logical Status
      Logical SCF, RIMP2
      External SCF, RIMP2
c$$$      Integer LUOut
c$$$      Parameter (LuOut = 6)
#include "stdio.fh"
*********************************************#include "stdio_init.fh"
      Status = SCF( RTDB)
      Status = RIMP2( RTDB)
      Return
      End
