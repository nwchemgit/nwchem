      SubRoutine ErrMem(LResult,Note,iRC)
*
* $Id: errmem.F,v 1.3 1997-11-04 10:10:38 d3e129 Exp $
*
      Implicit None
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
         Call ErrQuit(Note,iRC)
      End If
*
*---- Normal termination
      End
