      integer function sym_center_map(geom,cent,opnum)
C$Id$
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "global.fh"
c
      integer geom,cent,opnum
      integer addr
c
      if (.not. geom_check_handle(geom,'sym_center_map'))
     $     call errquit('sym_center_map: invalid geom', geom, GEOM_ERR)
      if (.not. geom_check_cent(geom, 'sym_center_map', cent))
     $     call errquit('sym_center_map: invalid center', cent,
     &       GEOM_ERR)
      if (opnum.lt.0 .or. opnum.gt.sym_num_ops(geom))
     $     call errquit('sym_center_map: invalid opnum', opnum,
     &       INPUT_ERR)
c
      addr = sym_center_map_index(geom) + opnum - 1
c
*      if (cent.eq.59 .or. cent.eq.50) then
*         write(6,*) ga_nodeid(), ' in sym_center_map', cent,
*     $        sym_center_map_index(geom),
*     $        addr,
*     $        sym_num_ops(geom)
*         call util_flush(6)
*      endif
c
      sym_center_map = int_mb(addr + (cent-1)*sym_num_ops(geom))

*      if (cent.eq.59 .or. cent.eq.50) then
*         write(6,*) ga_nodeid(), ' leaving sym_center_map', cent,
*     $        sym_center_map
*         call util_flush(6)
*      endif

*      write(6,1) opnum, cent, sym_center_map
* 1    format(' op ',i3,' * ',i3,' -> ',i3)
c
      end
