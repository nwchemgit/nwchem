      subroutine geom_print_rtdb_ecce(rtdb)
*
* $Id: geom_print_ecce.F,v 1.5 1999-09-22 19:00:16 d3g681 Exp $
*
#include "geom.fh"
      integer rtdb, geom
c
      if (.not. geom_create(geom, 'geometry'))  call errquit
     $     ('geom_print_rtdb_ecce: geom_create failed',0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry')) call errquit
     $     ('geom_print_rtdb_ecce: geom_rtdb_load failed',0)
      call geom_print_ecce(geom)
      if (.not. geom_destroy(geom)) call errquit
     $     ('geom_print_rtdb_ecce: geom destroy failed',0)
c
      end
      subroutine geom_print_ecce(geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
      integer geom
c
c     Output the tags, charges, and cartesian coordinates (Angs) 
c     to the ecce output file
c
      logical geom_check_handle
      external geom_check_handle
c
      integer nat
c
      if (.not. geom_check_handle(geom,'geom_print_ecce')) return
c
      nat = ncenter(geom)
      call dscal(nat*3, 1.0d0/angstrom_to_au, coords(1,1,geom), 1)
      call ecce_print2('cartesian coordinates', mt_dbl,
     $     coords(1,1,geom), 3, 3, nat)
      call dscal(nat*3, angstrom_to_au, coords(1,1,geom), 1)
c
      call ecce_print1('atomic charges', mt_dbl, charge(1,geom), nat)
      call ecce_print1_char('atomic tags', tags(1,geom), nat)
      call ecce_print1_char('group name', group_name(geom), 1)
c
      end

      
      
