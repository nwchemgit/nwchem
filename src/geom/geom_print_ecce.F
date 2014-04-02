      subroutine geom_print_rtdb_ecce(rtdb)
      implicit none
*
* $Id$
*
#include "geom.fh"
#include "errquit.fh"
      integer rtdb, geom
c
      if (.not. geom_create(geom, 'geometry'))  call errquit
     $     ('geom_print_rtdb_ecce: geom_create failed',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry')) call errquit
     $     ('geom_print_rtdb_ecce: geom_rtdb_load failed',0, RTDB_ERR)
      call geom_print_ecce(geom)
      if (.not. geom_destroy(geom)) call errquit
     $     ('geom_print_rtdb_ecce: geom destroy failed',0, GEOM_ERR)
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

*     **** output amatrix in angstroms ****
      if (isystype(geom).eq.3) then
         call dscal(3*3, 1.0d0/angstrom_to_au, amatrix(1,1,geom), 1)
         call ecce_print2('lattice vectors', mt_dbl,
     $     amatrix(1,1,geom), 3, 3, 3)
         call dscal(3*3, angstrom_to_au, amatrix(1,1,geom), 1)
      end if
c

      call ecce_print1('atomic charges', mt_dbl, charge(1,geom), nat)
      call ecce_print1_char('atomic tags', tags(1,geom), nat)
      call ecce_print1_char('group name', group_name(geom), 1)
c
      end

      
      
