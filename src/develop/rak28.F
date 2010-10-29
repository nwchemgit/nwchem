*
* $Id$
*
      logical function raktask_geomcalc(rtdb)
*test geom calc stuff
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
*::passed
      integer rtdb
*::local
      integer geom
*::code
      raktask_geomcalc = .false.
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911, RTDB_ERR)
      if (.not.geom_print_distances(geom)) stop ' dist die '
      if (.not.geom_print_angles(geom))    stop ' angl die '
      if (.not.geom_print_dihedrals(geom)) stop ' dihd die '
      raktask_geomcalc = geom_destroy(geom)
      end
