      SubRoutine ErrMem(LResult,Note,iRC)
*
* $Id$
*
      Implicit None
#include "errquit.fh"
*
************************************************************************
*
*---- Arguments of a ErrMem call
      Logical LResult
      Character*(*) Note
      Integer iRC
*
************************************************************************
*
      If (LResult) Then
         Return
      Else
         Call ErrQuit(Note,iRC, UNKNOWN_ERR)
      End If
*
*---- Normal termination
      End
