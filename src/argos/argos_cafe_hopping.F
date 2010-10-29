      logical function argos_cafe_hopping(lpbc,lpbcs,stimei,
     + isl,issgm,isgan,isq3,ishop,xs,nsaloc)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical argos_cafe_hop
      external argos_cafe_hop
c
      logical lpbc,lpbcs
      integer isl(msa,mis2)
      integer issgm(msa),nsaloc,isgan(msa),ishop(msa),isq3(msa)
      real*8 xs(msa,3)
      real*8 stimei
c
      integer i_itmp,l_itmp,i_dtmp,l_dtmp
c
      stime=stimei
      argos_cafe_hopping=.false.
c
      if(lpair.and.nhop.gt.0) then
      lhop=.true.
c
      if(.not.ma_push_get(mt_int,np,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',me)
c
      call argos_cafe_hoplist(issgm,nsaloc,int_mb(i_itmp),int_mb(i_lda),
     + dbl_mb(i_rda))
c
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to deallocate itmp',me)
c
      argos_cafe_hopping=.true.
      elseif(lhop) then
c
c     evaluate the E12
c
      call argos_cafe_uhop(lpbc,lpbcs,int_mb(i_lseq),isl,
     + issgm,isgan,isq3,xs,nsaloc,
     + dbl_mb(i_chg),int_mb(i_lda),dbl_mb(i_rda),dbl_mb(i_uda))
c
      if(.not.ma_push_get(mt_dbl,nldat,'dtmp',l_dtmp,i_dtmp))
     + call md_abort('Failed to allocate dtmp',me)
      if(.not.ma_push_get(mt_int,nldat,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',me)
c
      argos_cafe_hopping=argos_cafe_hop(int_mb(i_lseq),issgm,
     + int_mb(i_lda),dbl_mb(i_rda),dbl_mb(i_uda),
     + int_mb(i_itmp),dbl_mb(i_dtmp),isgan,ishop,xs,nsaloc,
     + int_mb(i_lsthop),dbl_mb(i_timhop))
c
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to deallocate itmp',me)
      if(.not.ma_pop_stack(l_dtmp))
     + call md_abort('Failed to deallocate dtmp',me)
c
      endif
c
      return
      end
