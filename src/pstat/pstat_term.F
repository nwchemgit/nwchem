      Logical Function PStat_Terminate()
C$Id: pstat_term.F,v 1.1 1995-12-13 01:36:36 d3g681 Exp $
      Implicit NONE
C
#include "pstatP.fh"
C
C     Nothing to do here, really, but let us "unassign" all of the items
C
      integer i
c
      do i = 1, pstat_max_items
         pstat_active(i) = .false.
      enddo
C
      PStat_Terminate = .TRUE.
      Return
      End
