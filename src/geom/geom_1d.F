      subroutine geom_1d(geom,c2au)
*
* $Id: geom_1d.F,v 1.3 1997-11-04 10:04:22 d3e129 Exp $
*
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "geomP.fh"
c
      integer geom
      double precision c2au
      write(*,*) 'SYSTYPE=POLYMER IS NOT A VALID OPTION AT THIS TIME'
c
      call errquit('geom_1d:invalid system type',isystype(1))
c
      return
      end










