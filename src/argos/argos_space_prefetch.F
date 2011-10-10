      subroutine argos_space_prefetch(n,iwl,xw,isl,xs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(2),ps(2)
      real*8 pwp(2),psp(2)
      logical lfetch,lnext
      integer ibbl,nfetched
c
      n=nbbl
c
      nfetched=0
      if(ibget.lt.0) nbget=nbbl+1
c
      do 1 ibbl=1,nbbl
      call argos_space_ibbl_prefetch(ibbl,int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl),iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),
     + xs,ps,psp,.false.,lfetch,lnext)
      icbbl=ibbl
      if(lfetch) nfetched=nfetched+1
      if(nfetched.ge.nbget.and..not.lnext) return
 1    continue
c
      return
      end
c $Id$
