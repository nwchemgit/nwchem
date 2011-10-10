      subroutine argos_space_ibbl_prefetch(ibbl,lbbl,ipl,jpl,
     + iwl,iwlp,xw,pw,pwp,isl,islp,xs,ps,psp,lpp,lfetch,lnext)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer ibbl
      integer lbbl(mbbl,mbb2),ipl(mbox,mip2),jpl(mbox,mip2)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp,lfetch,lnext
c
      integer i,j,indexw,indexs,indx,jndx,nndexw,nndexs
      integer jproc,jcell
      integer jccell,iccell
      integer ih,il,jh,jl
      integer nsnon,nwnon
      integer jwfr,jwto,jsfr,jsto
      integer handle_is,handle_iw,handle_s,handle_w
      integer handle_ps,handle_pw,handle_psp,handle_pwp
      character*80 string
c
      lfetch=.false.
c
      if(ibbl.eq.1) then
      jccell=-1
      else
      jccell=lbbl(ibbl-1,2)
      endif
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
c
      handle_iw=0
      handle_is=0
      handle_w=0
      handle_s=0
      handle_pw=0
      handle_ps=0
      handle_pwp=0
      handle_psp=0
c
      jproc=lbbl(ibbl,1)
      jcell=lbbl(ibbl,2)
c
      if(jproc.eq.me) then
c
      jccell=jcell
c
      else
c
      if(jcell.ne.jccell) then
c
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
c
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(indexw+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',indexw+nwnon)
      call ga_distribution(ga_iw,jproc,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,
     + iwl(indexw,1),mwm,handle_iw)
      else
      call ga_nbget(ga_iw,il+jwfr-1,il+jwto-1,jl,jl+npackw-1,
     + iwlp(indexw,1),mwm,handle_iw)
      endif
      call ga_distribution(ga_w,jproc,il,ih,jl,jh)
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(indexw,1,1),mwm,handle_w)
      if(lpp) then
      call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,jl+12*mwa+3,jl+18*mwa+2,
     + pw(indexw,1,1,1),mwm,handle_pw)
      if(lfree) call ga_nbget(ga_w,il+jwfr-1,il+jwto-1,
     + jl+18*mwa+3,jl+30*mwa+2,pwp(indexw,1,1,1,1),mwm,handle_pwp)
      endif
      lfetch=.true.
      endif
c
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(indexs+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (2)',indexs+nsnon)
      call ga_distribution(ga_is,jproc,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jh,
     + isl(indexs,1),msa,handle_is)
      else
      call ga_nbget(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(indexs,1),msa,handle_is)
      endif
      call ga_distribution(ga_s,jproc,il,ih,jl,jh)
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(indexs,1),msa,handle_s)
      if(lpp) then
      call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+12,jl+14,
     + ps(indexs,1,1),msa,handle_ps)
      if(lfree) call ga_nbget(ga_s,il+jsfr-1,il+jsto-1,jl+15,jl+20,
     + psp(indexs,1,1,1),msa,handle_psp)
      endif
      lfetch=.true.
      endif
c
      endif
c
      jccell=jcell
c
      endif
c
      lbbl(ibbl,15)=handle_iw
      lbbl(ibbl,16)=handle_w
      lbbl(ibbl,17)=handle_is
      lbbl(ibbl,18)=handle_s
      lbbl(ibbl,19)=handle_pw
      lbbl(ibbl,20)=handle_ps
      lbbl(ibbl,21)=handle_pwp
      lbbl(ibbl,22)=handle_psp
c
      lnext=.false.
      if(ibbl.lt.nbbl) then
      lnext=lbbl(ibbl+1,2).eq.lbbl(ibbl,2)
      endif
c
      return
      end
c $Id$
