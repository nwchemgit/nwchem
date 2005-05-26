      subroutine ffflush(unit)
C$Id: ffflush.F,v 1.4 2005-05-26 00:23:57 edo Exp $
      integer unit              ! [input]
c
c     Flush outputto this fortran unit
c
      call util_flush(unit)
c
      end
