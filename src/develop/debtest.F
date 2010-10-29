      Logical Function DEBTest(RTDB)
C$Id$
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
