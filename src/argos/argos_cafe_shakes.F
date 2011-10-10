      subroutine argos_cafe_shakes(lhandl,xs,ys,igan,isat,
     + isgm,isdt,ishp,isfr,isto)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lhandl,isfr,isto
      real*8 xs(msa,3),ys(msa,3)
      integer igan(msa),isat(msa),isgm(msa),isdt(msa),ishp(msa)
c
      integer lptr,nbonds
c
      lptr=int_mb(i_list+lhandl+4)
      nbonds=int_mb(lptr)
      if(nbonds.gt.0) then
      call argos_cafe_shks(nbonds,int_mb(lptr+1),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + isto,msa,igan,isat,isgm,isdt,ishp,xs,ys,dbl_mb(i_mas),
     + int_mb(i_lseq))
      endif
c
      return
      end
c $Id$
