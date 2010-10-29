      subroutine geom_1d(geom,c2au)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      double precision c2au
      write(*,*) 'SYSTYPE=POLYMER IS NOT A VALID OPTION AT THIS TIME'
c
      call errquit('geom_1d:invalid system type',0,GEOM_ERR)
c
      return
      end










