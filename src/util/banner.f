      Subroutine Banner(LU, Msg, Char, Top, Bot, Sides)
C$Id$
      Implicit NONE
      Integer LU
      Character*(*) Msg
      Character*(1) Char
      Logical Top, Bot, Sides
C
      Character*(80) Fmt
      Integer MsgLen
C
C     Figure out message length & create an appropriate format
C
      MsgLen = Len(Msg)
C     
      If ( Sides ) then
         Write( Fmt, 9000) MsgLen+4, Char
      Else
         Write( Fmt, 9000) MsgLen, Char
      EndIf
 9000 Format('(1X,', I5, '(''', A, '''))')
C
      If ( Top ) Write (Lu, Fmt=Fmt)
      If ( Sides ) then
         Write (Lu, Fmt = 9010) Char, Msg, Char
      Else
         Write (Lu, Fmt = 9010) Msg
      EndIf
      If ( Bot ) Write (Lu, Fmt=Fmt)
 9010 Format(1X, A, 1X, A, 1X, A)
C
      Return
      End
