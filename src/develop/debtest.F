      Subroutine DEBTest(RTDB)
C$Id: debtest.F,v 1.3 1995-02-02 23:12:34 d3g681 Exp $
      Implicit NONE
      Integer RTDB
      Logical Status
      Logical SCF, RIMP2
      External SCF, RIMP2
c$$$      Integer LUOut
c$$$      Parameter (LuOut = 6)
#include "stdio.fh"
#include "stdio_init.fh"
      Status = SCF( RTDB)
      Status = RIMP2( RTDB)
      Return
      End
