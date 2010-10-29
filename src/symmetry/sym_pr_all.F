      subroutine sym_print_all(geom, oinfo, ouniq, omap, oops, ochar)
C$Id$
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "inp.fh"
c
c     Print out all symmetry related information inside the
c     geometry object 
c
      integer geom              ! [input]
      logical oinfo             ! If true print info
      logical ouniq             ! If true print list of unique atoms
      logical omap              ! If true print mapping of atoms under ops
      logical oops              ! If true print operator matrices
      logical ochar             ! If true print character table
c
      integer cent, j, op, addr, centhi
c
      call util_print_centered(6,'Symmetry information', 15, .true.)
      write(6,*)
c
      if (oinfo) then
         write(6,1) group_name(geom), group_number(geom),
     $        sym_num_ops(geom)+1, ncenter_unique(geom)
         if (isystype(geom) .ne. 0) write(6,11) 
     $        setting_number(geom)
         write(6,*)
 1       format(
     $        ' Group name             ', a10/
     $        ' Group number           ', i4/
     $        ' Group order            ', i4/
     $        ' No. of unique centers  ', i4)
 11      format(
     $        ' Setting number         ', i4)
      endif
c      
      if (ochar) call sym_print_char_table(geom)
c
      if (ouniq) then
         call util_print_centered(6,'Symmetry unique atoms', 15,.false.)
         write(6,*)
         do j = 1, ncenter_unique(geom), 15
            centhi = min(j + 14,ncenter_unique(geom))
            write(6,2) (unique_cent(cent,geom),cent=j,centhi)
 2          format(1x,15i5)
         enddo
         write(6,*)
      endif
c
      if (omap) then
         call util_print_centered(6,'Mapping of centers', 14, .false.)
         write(6,*)
         do op = 1, sym_num_ops(geom)
            write(6,15) op
 15         format('     op ', i3)
            if (omap) then
               addr = sym_center_map_index(geom) + op - 1
               write(6,5) (
     $              cent,
     $              int_mb(addr + (cent-1)*sym_num_ops(geom)),
     $              cent=1,ncenter(geom))
 5             format(100(10x,5('(',i4,'->',i4,')':,2x)/))
            endif
         enddo
         write(6,*)
      endif
c
      if (oops) then
         call util_print_centered(6,'Group operations', 13, .false.)
         write(6,*)
         call sym_print_ops(geom)
         write(6,*)
      endif
c
      call util_flush(6)
c
      end
