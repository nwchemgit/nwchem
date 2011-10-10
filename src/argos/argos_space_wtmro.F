      subroutine argos_space_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,iwlp,xw,vw,xwcr,isl,islp,xs,vs,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer ipl(mbox,mip2)
c
      integer j,k,l,nwmn,nsan,node,ilp,ihp,jlp,jhp,ili,ihi,jli,jhi
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      character*10 rdate,rtime
      character*18 user
#ifdef USE_POSIXF
      integer ilen,ierror
#endif
c
      write(lfnmro) nwm,nwa,nsa,stime,temp,pres,vlat,nhist
c
      call swatch(rdate,rtime)
#ifdef USE_POSIXF
      call pxfgetlogin(user, ilen, ierror)
#else
      call getlog(user)
#endif
      if(user(18:18).ne.' ') user='                  '
      hist(nhist)(1:18)=user
      hist(nhist)(19:28)=rdate
      hist(nhist)(29:48)=rtime
      hist(nhist)(49:108)=project(1:60)
      write(lfnmro) (hist(j),j=1,nhist)
c
      do 1 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      write(lfnmro) ((ipl(j,k),j=1,mbox),k=1,mip2)
      nwmn=ipl(1,2)
      nsan=ipl(2,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jhi,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call argos_space_unpackw(nwmn,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      write(lfnmro) ((iwl(j,k),j=1,nwmn),k=1,miw2)
      write(lfnmro) (((xw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      write(lfnmro) (((vw(j,k,l),j=1,nwmn),k=1,3),l=1,nwa)
      endif
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jhi,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call argos_space_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      write(lfnmro) ((isl(j,k),j=1,nsan),k=1,mis2)
      write(lfnmro) ((xs(j,k),j=1,nsan),k=1,3)
      write(lfnmro) ((vs(j,k),j=1,nsan),k=1,3)
      endif
    1 continue
c
      return
      end
c $Id$
