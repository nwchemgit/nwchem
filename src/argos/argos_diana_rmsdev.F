      subroutine argos_diana_rmsdev(xtrj,xr,wt,wts,xrms,rms,trms,ga_rms,
     + isel,sgmnam)
c $Id$
      implicit none
c
#include "argos_diana_params.fh"
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ga_rms
      real*8 xtrj(ntrj,3,nfrdim)
      real*8 xr(nsa,3),wt(nsa),wts(nsa),xrms(nsa)
      real*8 rms(2,nfrdim),trms(3,nfrtot)
      character*16 sgmnam(nsa)
      integer isel(nsa)
c
      real*8 arms,bfact
c
      integer i,j,k,n,ip,ilc,ihc,jlc,jhc,ilt,iht,jlt,jht
c
      bfact=128.0d0*atan(1.0d0)*atan(1.0d0)/3.0d0
c
      do 1 k=1,nfrme
      call super(xtrj(1,1,k),ntrj,ntrj,xr,wt,wts,xrms,nsa,nsa,
     + .false.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.false.)
      rms(1,k)=rms0
      rms(2,k)=rms1
    1 continue
c
      call ga_dgop(mag_d01,xrms,nsa,'+')
c
      if(me.eq.0) then
c
      open(unit=lfnrms,file=filrms(1:index(filrms,' ')-1),
     + form='formatted',status='unknown')
c
      nsgm=1
      do 2 i=2,nsa
      if(sgmnam(i)(11:16).ne.sgmnam(i-1)(11:16)) nsgm=nsgm+1
    2 continue
c
      write(lfnrms,1000) nsa,nsgm,nfrtot
 1000 format(3i7)
c
      if(me.eq.0) then
c
      do 3 i=1,nsa
      xrms(i)=xrms(i)/dble(nfrtot)
    3 continue
c
      write(lfnrms,1001) (sgmnam(i),i,isel(i),
     + sqrt(xrms(i)),bfact*xrms(i),i=1,nsa)
 1001 format(a16,i10,i5,2f12.6)
c
      n=0
      do 4 i=1,nsa
      if(isel(i).ne.0) then
      arms=arms+xrms(i)
      n=n+1
      endif
      if(i.eq.nsa) goto 5
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 5
      goto 4
    5 continue
      if(n.gt.0) then
      write(lfnrms,1002) sgmnam(i)(1:5),sgmnam(i)(11:16),
     + sqrt(arms/dble(n)),bfact*arms/dble(n)
 1002 format(a5,a6,2f12.6)
      endif
      arms=0.0d0
      n=0
    4 continue
c
      endif
c
      do 6 ip=0,np-1
      call ga_distribution(ga_rms,ip,ilc,ihc,jlc,jhc)
      call ga_get(ga_rms,1,2,jlc,jhc,trms(2,ip*nfrdim+1),3)
      call ga_distribution(ga_trt,ip,ilt,iht,jlt,jht)
      call ga_get(ga_trt,1,1,jlt,jht,trms(1,ip*nfrdim+1),3)
    6 continue
c
      write(lfnrms,1003) ((trms(j,k),j=1,3),k=1,nfrtot)
 1003 format(3f12.6)
c
      close(unit=lfnrms,status='keep')
      endif
c
      return
      end
