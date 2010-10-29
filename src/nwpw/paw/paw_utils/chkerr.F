c
c $Id$
c

      subroutine chkerr(string, icode)
      implicit none

      character*(*) string
      integer icode
      
      if(icode.ne.0) then
        print*,"aborting due to error ",string
      stop
      end if

      end
