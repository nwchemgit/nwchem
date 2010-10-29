      subroutine ffflush(unit)
C$Id$
      integer unit              ! [input]
c
c     Flush outputto this fortran unit
c
      call util_flush(unit)
c
      end
