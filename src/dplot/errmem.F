      SubRoutine ErrMem(LResult,Note,iRC)
*
* $Id: errmem.F,v 1.4 2003-10-17 22:54:43 carlfahl Exp $
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
