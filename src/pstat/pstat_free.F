      Logical Function PStat_Free( Handle )
      Implicit NONE
      Integer Handle
C
#include "pstatP.fh"
C
C     For now just blank out the name
C
      PStat_Name( Handle ) = 'Freed!'
C
      PStat_Free = .TRUE.
      Return
      End
