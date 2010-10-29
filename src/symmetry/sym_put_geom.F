      subroutine sym_put_in_geom(geom, nat_uniq, nat_new,
     $     tags_new, uniq_cent, coords_new, l_map, k_map, nops,
     $     charge_new, mass_new, invnucexp_new)
C$Id$
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "geomP.fh"

      integer geom,nat_uniq,nat_new,i,l_map, k_map, nops
      character*16 tags_new(nat_new)
      integer uniq_cent(nat_uniq)
      double precision coords_new(3,nat_new)
      double precision charge_new(max_cent)
      double precision mass_new(max_cent)
      double precision invnucexp_new(max_cent)
c
c     tags_new(1:nat_new) = tags for full list of centers
c     nat_uniq = no. of symmetry unique atoms = no. of atoms in the geom
c     nat_new = no. of atoms in full list
c     coords_new = list of FRACTIONAL coordinates that must be
c                  converted to CARTESIAN coordinates before insertion
c     l_map/k_map = MA handle/index to atom symmetry map
c     nop = no. of operations in the group
c
      do i = 1, nat_new
         charge(i,geom) = charge_new(i)
         geom_mass(i,geom) = mass_new(i)
         geom_invnucexp(i,geom) = invnucexp_new(i)
         tags(i,geom)=tags_new(i)
      enddo
      call dgemm('n', 'n', 3, nat_new, 3, 1.0d0, amatrix(1,1,geom), 3,
     $     coords_new, 3, 0.0d0, coords(1,1,geom), 3)
c
      ncenter_unique(geom) = nat_uniq
      ncenter(geom) = nat_new
      do i=1,nat_uniq
        unique_cent(i,geom)=uniq_cent(i)
      enddo
c
      if (sym_center_map_handle(geom) .ne. -1) then
         if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $        call errquit('sym_nwc: free of old atom map', 0, MA_ERR)
      endif
      sym_center_map_handle(geom) = l_map
      sym_center_map_index(geom)  = k_map
      sym_num_ops(geom) = nops
c
      end
