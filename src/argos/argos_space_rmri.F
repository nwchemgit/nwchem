      logical function argos_space_rmri(lfnmri,stime,pres,
     + temp,tempw,temps,
     + iwl,iwlp,xw,vw,xwcr,isl,islp,xs,vs,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfnmri
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer i,j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ltemp(3)
      real*8 rtemp(12)
c
      if(me.eq.0) then
      read(lfnmri,err=9,end=9) ltemp,rtemp,nhist
      if(nhist.gt.0) read(lfnmri,err=9,end=9) (hist(j),j=1,nhist)
      do 1 node=np-1,0,-1
      read(lfnmri,err=9) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_put(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      if(nwmn.gt.0) then
      read(lfnmri,err=9) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      read(lfnmri,err=9) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      read(lfnmri,err=9) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_put(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      else
      call argos_space_packw(nwmn,iwl,iwlp)
      call ga_put(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsan.gt.0) then
      read(lfnmri,err=9) ((isl(j,k),j=1,nsan),k=1,mis2)
      read(lfnmri,err=9) ((xs(j,k),j=1,nsan),k=1,3)
      read(lfnmri,err=9) ((vs(j,k),j=1,nsan),k=1,3)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_put(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      else
      call argos_space_pack(nsan,isl,islp)
      call ga_put(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_put(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_put(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      endif
    1 continue
      endif
c
      call ga_brdcst(msp_12,ltemp,3*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(msp_13,rtemp,12*ma_sizeof(mt_dbl,1,mt_byte),0)
      nwm=ltemp(1)
      nwa=ltemp(2)
      nsa=ltemp(3)
      stime=rtemp(1)
      temp=rtemp(2)
      pres=rtemp(3)
      k=3
      do 2 i=1,3
      do 3 j=1,3
      k=k+1
      vlat(i,j)=rtemp(k)
    3 continue
    2 continue
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmloc=ipl(1,2)
      nsaloc=ipl(2,2)
      if(nwmloc.gt.0) then
      call ga_distribution(ga_iw,me,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmloc-1,jli,jhi,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmloc-1,jli,jli+npackw-1,iwlp,mwm)
      call argos_space_unpackw(nwmloc,iwl,iwlp)
      endif
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw,jlw+3*mwa-1,xw,mwm)
      call ga_get(ga_w,ilw,ilw+nwmloc-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      endif
      if(nsaloc.gt.0) then
      call ga_distribution(ga_is,me,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsaloc-1,jli,jhi,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsaloc-1,jli,jli+npack-1,islp,msa)
      call argos_space_unpack(nsaloc,isl,islp)
      endif
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsaloc-1,jls+3,jls+5,vs,msa)
      endif
c
      argos_space_rmri=.true.
      return
c
    9 continue
      argos_space_rmri=.false.
      return
      end
c $Id$
