      subroutine argos_space_printf(filtop,lfntop,isl,xs,fs,
     + npener,esa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfntop
      character*70 filtop
      integer isl(msa,mis2),npener
      real*8 xs(msa,3),fs(msa,3,2),esa(nsa,2)
c
      integer lenscr
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((6*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call md_abort('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call md_abort('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call md_abort('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call md_abort('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call md_abort('Failed to allocate n',0)
c
      call argos_space_prt_s(filtop,lfntop,
     + int_mb(i_ipl),isl,int_mb(i_pack),xs,fs,
     + lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + npener,esa)
c
      if(.not.ma_pop_stack(l_n))
     + call md_abort('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call md_abort('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call md_abort('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bv))
     + call md_abort('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call md_abort('Failed to deallocate bx',0)
c
      return
      end
c $Id$
