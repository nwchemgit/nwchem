c
c $Id: chg_case.f,v 1.1 2008-04-18 17:40:29 marat Exp $
c
      SUBROUTINE chg_case(char)

      implicit none

      integer upper_to_lower

      character*1 char

      upper_to_lower=iachar("a")-iachar("A")

      if("a"<=char.and.char<="z")then
        char=ACHAR(IACHAR(char)-upper_to_lower)
      endif

      return

      END
