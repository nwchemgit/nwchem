c
c $Id: chkerr.F,v 1.1 2003-04-07 21:52:36 windus Exp $
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
