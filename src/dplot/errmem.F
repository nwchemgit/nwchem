      SubRoutine ErrMem(LResult,Note,iRC)
      Implicit None
*
************************************************************************
*
*---- Arguments of a ErrMem call
      Logical LResult
      Character*80 Note
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
