      subroutine argos_space_alloc
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer isize
c
      if(.not.ma_push_get(mt_int,np,'bindex',l_bindex,i_bindex))
     + call md_abort('Failed to allocate bindex',0)
      if(.not.ma_push_get(mt_int,54*np,'buren',l_buren,i_buren))
     + call md_abort('Failed to allocate buren',0)
      if(.not.ma_push_get(mt_int,3*maxbox,'owner',l_iown,i_iown))
     + call md_abort('Failed to allocate owner',0)
      if(.not.ma_push_get(mt_dbl,3*maxbox,'bxsiz',l_boxs,i_boxs))
     + call md_abort('Failed to allocate bxsiz',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'ibxrg',l_boxr,i_boxr))
     + call md_abort('Failed to allocate ibxrg',0)
      if(.not.ma_push_get(mt_int,6*maxbox,'rng',l_rng,i_rng))
     + call md_abort('Failed to allocate rng',me)
c
      if(.not.ma_push_get(mt_int,mip2*mbox,'ipl',l_ipl,i_ipl))
     + call md_abort('Failed to allocate ipl',0)
      if(.not.ma_push_get(mt_int,mip2*mbox,'jpl',l_jpl,i_jpl))
     + call md_abort('Failed to allocate jpl',0)
c
      mseq=nseq
c
      if(.not.ma_push_get(mt_int,mseq,'lseq',l_lseq,i_lseq))
     + call md_abort('Failed to allocate lseq',0)
      if(.not.ma_push_get(mt_int,mseq,'sndx',l_sndx,i_sndx))
     + call md_abort('Failed to allocate sndx',0)
c
      if(.not.ga_create(mt_int,np*mbox,mip2,'ip',mbox,mip2,ga_ip))
     + call md_abort('Failed to create global array ip',0)
c
      return
      end
c $Id$
