      subroutine sym_group_name(geom, name)
C$Id: sym_grp_name.F,v 1.2 1999-07-14 00:29:11 d3e129 Exp $
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom              ! [input]
      character*(*) name        ! [output]
c
c     return the name of the group
c
      integer i
      logical geom_check_handle
      external geom_check_handle
c
      if (.not. geom_check_handle(geom,'sym_group_name'))
     $     call errquit('sym_group_name: bad geom handle',geom)
c
c     strip annonying leading blanks
c
      do i = 1, len(group_name(geom))
         if (group_name(geom)(i:i).ne.' ') goto 10
      enddo
      call errquit('sym_group_name: all blank?', 0)
c
 10   name = group_name(geom)(i:)
c
      end
