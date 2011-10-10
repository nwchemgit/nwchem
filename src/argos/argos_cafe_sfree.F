      subroutine argos_cafe_sfree()
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
c      if(.not.lscr) call md_abort('Error 1 in argos_cafe_sfree',0)
      if(.not.lscr) return
c
      if(ipme.gt.0) then
      if(.not.ma_free_heap(l_dtheta))
     + call md_abort('Failed to free scratch array dtheta',0)
      if(.not.ma_free_heap(l_theta))
     + call md_abort('Failed to free scratch array theta',0)
      if(.not.ma_free_heap(l_pmew))
     + call md_abort('Failed to free scratch array pmew',0)
      if(.not.ma_free_heap(l_pmes))
     + call md_abort('Failed to free scratch array pmes',0)
      endif
c
      if(.not.ma_free_heap(l_smr6))
     + call md_abort('Failed to free scratch array smr6',0)
c
      if(.not.ma_free_heap(l_smr5))
     + call md_abort('Failed to free scratch array smr5',0)
c
      if(.not.ma_free_heap(l_smr4))
     + call md_abort('Failed to free scratch array smr4',0)
c
      if(.not.ma_free_heap(l_smr3))
     + call md_abort('Failed to free scratch array smr3',0)
c
      if(.not.ma_free_heap(l_smr2))
     + call md_abort('Failed to free scratch array smr2',0)
c
      if(.not.ma_free_heap(l_smr1))
     + call md_abort('Failed to free scratch array smr1',0)
c
      if(.not.ma_free_heap(l_s3r2))
     + call md_abort('Failed to free scratch array s3r2',0)
c
      if(.not.ma_free_heap(l_s3r1))
     + call md_abort('Failed to free scratch array s3r1',0)
c
      if(.not.ma_free_heap(l_s1r6))
     + call md_abort('Failed to free scratch array s1r6',0)
c
      if(.not.ma_free_heap(l_s1r5))
     + call md_abort('Failed to free scratch array s1r5',0)
c
      if(.not.ma_free_heap(l_s1r4))
     + call md_abort('Failed to free scratch array s1r4',0)
c
      if(.not.ma_free_heap(l_s1r3))
     + call md_abort('Failed to free scratch array s1r3',0)
c
      if(.not.ma_free_heap(l_s1r2))
     + call md_abort('Failed to free scratch array s1r2',0)
c
      if(.not.ma_free_heap(l_s1r1))
     + call md_abort('Failed to free scratch array s1r1',0)
c
      if(.not.ma_free_heap(l_s2i1))
     + call md_abort('Failed to free scratch array s2i1',0)
c
      if(.not.ma_free_heap(l_s1i6))
     + call md_abort('Failed to free scratch array s1i6',0)
c
      if(.not.ma_free_heap(l_s1i5))
     + call md_abort('Failed to free scratch array s1i5',0)
c
      if(.not.ma_free_heap(l_s1i4))
     + call md_abort('Failed to free scratch array s1i4',0)
c
      if(.not.ma_free_heap(l_s1i3))
     + call md_abort('Failed to free scratch array s1i3',0)
c
      if(.not.ma_free_heap(l_s1i2))
     + call md_abort('Failed to free scratch array s1i2',0)
c
      if(.not.ma_free_heap(l_s1i1))
     + call md_abort('Failed to free scratch array s1i1',0)
c
      lscr=.false.
c
      return
      end
c $Id$
