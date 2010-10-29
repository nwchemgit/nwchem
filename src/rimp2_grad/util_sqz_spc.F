C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     util_sqz_spc -- Squeeze multiple spaces from strings
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Util_Sqz_Spc(String)
      Implicit NONE
C
C ARGUMENTS
      Character*(*) String ![inout] String to squeeze
C
C DESCRIPTION
C     Squeezes multiple spaces to a single space.  Eliminates ALL
C     leading blanks.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Src, Dest
      Logical First
C
C     Starting with First true eliminates all leading blanks.  False
C     will leave one.
C
      First = .TRUE.
      Dest = 1
      Do Src = 1, Len( String )
         If ( .NOT. First .OR. String(Src:Src) .ne. ' ') then
            String(Dest:Dest) = String(Src:Src)
            Dest = Dest + 1
            If ( String(Src:Src) .eq. ' ') then
               First = .TRUE.
            Else
               First = .FALSE.
            EndIf
         EndIf
      EndDo
C
C     Now make sure the rest is filled with blanks.  Dest already
C     advanced to next position by loop above.
C
      Do Src = Dest, Len( String )
         String(Src:Src) = ' '
      EndDo
C
      Return
      End
c$$$      Program Test
c$$$      Character*40 Input
c$$$C
c$$$ 1000 Continue
c$$$C
c$$$      Write (6, *) 'Enter string:'
c$$$      Read (5, '(A)') Input
c$$$      Write (6, 9000) 'Input', Input
c$$$      Call Util_Sqz_Spc(Input)
c$$$      Write (6, 9000) 'Output', Input
c$$$ 9000 Format(A6, 2X, '|', A, '|')
c$$$C
c$$$      Goto 1000
c$$$C
c$$$      Stop
c$$$      End
