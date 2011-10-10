      subroutine argos_space_gtxpbl(ibbl,lhandl,
     + iwl,iwlp,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,islp,xs,ps,psp,isfr,isto,jsfr,jsto,
     + lbbl,ipl,jpl,lself,lpp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2),jpl(mbox,mip2)
      logical lself,lpp
c
      integer i,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      character*80 string
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
c
      lself=ibox.eq.jbox
c
      if(ibbl.eq.1) then
      jcnode=-1
      jcbox=-1
      icbox=-1
      endif
c
      if(ibox.ne.icbox) then
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
 1000 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1001) ibox,me
 1001 format('SP0001: Could not find local ibox',i5,' on node',i5)
      call md_abort(string,me)
      endif
c
      icbox=ibox
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
      endif
c
      if(jnode.eq.me) then
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 2 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.jbox) ndx=30+i
    2 continue
c
      if(ndx.eq.0) then
      write(*,1000) me,(ipl(30+i,1),i=1,ipl(1,1))
      write(string,1002) jbox,me
 1002 format('SP0001: Could not find local jbox',i5,' on node',i5)
      call md_abort(string,me)
      endif
c
      jcbox=jbox
      jwfr=ipl(ndx,2)
      jwto=ipl(ndx,3)
      jsfr=ipl(ndx,4)
      jsto=ipl(ndx,5)
      endif
c
      else
c
      if(jnode.ne.jcnode) then
      call ga_distribution(ga_ip,jnode,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,jpl,mbox)
      jcnode=jnode
      endif
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0) then
      if(iand(idebug,4).eq.4) then
      write(lfndbg,8000) jbox,jnode
 8000 format('Could not find box',i5,' on node',i5,
     + ' in argos_space_gtxpbl')
      write(lfndbg,8001) jnode,(jpl(30+i,1),i=1,jpl(1,1))
 8001 format('Box list obtained from node',i5,' is',/,(20i5))
      call util_flush(lfndbg)
      endif
c
      write(*,1003) jnode,(jpl(30+i,1),i=1,jpl(1,1))
 1003 format(/,'Boxes on node',i5,':',t20,20i5,/,(t20,20i5))
      write(string,1004) jbox,jnode
 1004 format('SP0002: Could not find remote jbox',i5,' on node',i5)
      call md_abort(string,me)
c
      endif
      jcbox=jbox
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_iw,jnode,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_get(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,iwl(nwmloc+1,1),mwm)
      else
      call ga_get(ga_iw,il+jwfr-1,il+jwto-1,jl,jl+npackw-1,
     + iwlp(nwmloc+1,1),mwm)
      call argos_space_unpackw(jwto-jwfr+1,iwl(nwmloc+1,1),
     + iwlp(nwmloc+1,1))
      endif
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_get(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(nwmloc+1,1,1),mwm)
      if(lpp) then
      call ga_get(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+18*mwa+2,
     + pw(nwmloc+1,1,1,1),mwm)
      if(lfree) call ga_get(ga_w,il+jwfr-1,il+jwto-1,
     + jl+18*mwa+3,jl+30*mwa+2,pwp(nwmloc+1,1,1,1,1),mwm)
      endif
      jwfr=nwmloc+1
      jwto=nwmloc+1+jpl(ndx,3)-jpl(ndx,2)
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (3)',nsaloc+nsnon)
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jh,isl(nsaloc+1,1),msa)
      else
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(nsaloc+1,1),msa)
      call argos_space_unpack(jsto-jsfr+1,isl(nsaloc+1,1),
     + islp(nsaloc+1,1))
      endif
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(nsaloc+1,1),msa)
      if(lpp) then
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(nsaloc+1,1,1),msa)
      if(lfree) call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(nsaloc+1,1,1,1),msa)
      endif
      jsfr=nsaloc+1
      jsto=nsaloc+1+jpl(ndx,5)-jpl(ndx,4)
      endif
      endif
c
      endif
c
      return
      end
c $Id$
