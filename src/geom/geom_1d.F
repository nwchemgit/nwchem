      subroutine geom_1d(geom)
c
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "geomP.fh"
c
      integer geom
      write(*,*) 'SYSTYPE=POLYMER IS NOT A VALID OPTION AT THIS TIME'
c
      call errquit('geom_1d:invalid system type',isystype(1))
c
      return
      end










