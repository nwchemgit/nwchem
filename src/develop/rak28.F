*
* $Id: rak28.F,v 1.2 1999-07-28 00:42:15 d3e129 Exp $
*
      logical function raktask_geomcalc(rtdb)
*test geom calc stuff
      implicit none
#include "stdio.fh"
#include "geom.fh"
*::passed
      integer rtdb
*::local
      integer geom
*::code
      raktask_geomcalc = .false.
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911)
      if (.not.geom_print_distances(geom)) stop ' dist die '
      if (.not.geom_print_angles(geom))    stop ' angl die '
      if (.not.geom_print_dihedrals(geom)) stop ' dihd die '
      raktask_geomcalc = geom_destroy(geom)
      end
