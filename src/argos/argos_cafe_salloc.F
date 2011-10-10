      subroutine argos_cafe_salloc()
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lenscr
c
      if(lscr) call md_abort('Error 1 in argos_cafe_salloc',0)
      if(mscr.le.0) call md_abort('Error 2 in argos_cafe_salloc',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i1',l_s1i1,i_s1i1))
     + call md_abort('Failed to allocate scratch array s1i1',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i2',l_s1i2,i_s1i2))
     + call md_abort('Failed to allocate scratch array s1i2',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i3',l_s1i3,i_s1i3))
     + call md_abort('Failed to allocate scratch array s1i3',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i4',l_s1i4,i_s1i4))
     + call md_abort('Failed to allocate scratch array s1i4',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i5',l_s1i5,i_s1i5))
     + call md_abort('Failed to allocate scratch array s1i5',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i6',l_s1i6,i_s1i6))
     + call md_abort('Failed to allocate scratch array s1i6',0)
c
      lenscr=max(nsatot,3*mscr)
      if(.not.ma_alloc_get(mt_int,lenscr,'s2i1',l_s2i1,i_s2i1))
     + call md_abort('Failed to allocate scratch array s2i1',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r1',l_s1r1,i_s1r1))
     + call md_abort('Failed to allocate scratch array s1r1',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r2',l_s1r2,i_s1r2))
     + call md_abort('Failed to allocate scratch array s1r2',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r3',l_s1r3,i_s1r3))
     + call md_abort('Failed to allocate scratch array s1r3',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r4',l_s1r4,i_s1r4))
     + call md_abort('Failed to allocate scratch array s1r4',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r5',l_s1r5,i_s1r5))
     + call md_abort('Failed to allocate scratch array s1r5',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r6',l_s1r6,i_s1r6))
     + call md_abort('Failed to allocate scratch array s1r6',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mscr,'s3r1',l_s3r1,i_s3r1))
     + call md_abort('Failed to allocate scratch array s3r1',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mscr,'s3r2',l_s3r2,i_s3r2))
     + call md_abort('Failed to allocate scratch array s3r2',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr1',l_smr1,i_smr1))
     + call md_abort('Failed to allocate scratch array smr1',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr2',l_smr2,i_smr2))
     + call md_abort('Failed to allocate scratch array smr2',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr3',l_smr3,i_smr3))
     + call md_abort('Failed to allocate scratch array smr3',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr4',l_smr4,i_smr4))
     + call md_abort('Failed to allocate scratch array smr4',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr5',l_smr5,i_smr5))
     + call md_abort('Failed to allocate scratch array smr5',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr6',l_smr6,i_smr6))
     + call md_abort('Failed to allocate scratch array smr6',0)
c
      if(ipme.gt.0) then
      if(.not.ma_alloc_get(mt_dbl,3*msa,'pmes',l_pmes,i_pmes))
     + call md_abort('Failed to allocate scratch array pmes',0)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'pmew',l_pmew,i_pmew))
     + call md_abort('Failed to allocate scratch array pmew',0)
      lenscr=3*morder*(mwa*mwm+msa)
      if(.not.ma_alloc_get(mt_dbl,lenscr,'theta',l_theta,i_theta))
     + call md_abort('Failed to allocate scratch array theta',0)
      if(.not.ma_alloc_get(mt_dbl,lenscr,'dtheta',l_dtheta,i_dtheta))
     + call md_abort('Failed to allocate scratch array dtheta',0)
      endif
c
      lscr=.true.
c
      return
      end
c $Id$
