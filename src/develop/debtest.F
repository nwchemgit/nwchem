      Logical Function DEBTest(RTDB)
C$Id: debtest.F,v 1.5 1996-03-15 22:04:26 gg502 Exp $
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
      DEBTest = Status
      Return
      End
