      subroutine ffflush(unit)
      integer unit              ! [input]
c
c     Flush outputto this fortran unit
c
      call flush(unit)
c
      end
