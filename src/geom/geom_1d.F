      subroutine geom_1d(geom,c2au)
*
* $Id: geom_1d.F,v 1.5 2003-10-17 22:55:06 carlfahl Exp $
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
      call errquit('geom_1d:invalid system type',isystype(1), GEOM_ERR)
c
      return
      end










