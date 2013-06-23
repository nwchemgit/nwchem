c
c     qmd_geom_store: store geometry
c
      subroutine qmd_geom_store(rtdb,geom,nat,r)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "geom.fh"
c
      integer rtdb               ! rtdb
      integer geom               ! geom
      integer nat                ! number of atoms
      double precision r(3,nat)  ! coordinates
c
      character*30 pname
c
c     preliminaries
      pname = "qmd_geom_store: "
c
c     update geom with coords for energy-gradient calculation
      if (.not.geom_cart_coords_set(geom,r))
     $   call errquit(pname//'geom_cart_coords_set',0,GEOM_ERR)
c
c     store updated geometry
      if(.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &   call errquit(pname//'geom_rtdb_store',0,GEOM_ERR)
c
      return
      end
c $Id$
