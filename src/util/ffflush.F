      subroutine ffflush(unit)
C$Id: ffflush.F,v 1.3 1995-02-02 18:09:23 d3g681 Exp $
      integer unit              ! [input]
c
c     Flush outputto this fortran unit
c
      call flush(unit)
c
      end
