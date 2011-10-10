      subroutine argos_space_nbwait2(ibbl,lnew,lbbl,lhandl,lself,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,iwlp,isl,islp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,lhandl
      integer lbbl(mbbl,mbb2)
      integer isl(msa,mis2),islp(msa,npack)
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      logical lself,lnew
      integer indexs,indexw
      integer handle_iw,handle_w,handle_is,handle_s
      integer handle_pw,handle_ps,handle_pwp,handle_psp
      integer i
c
      lnew=.false.
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
      iwfr=lbbl(ibbl, 7)
      iwto=lbbl(ibbl, 8)
      isfr=lbbl(ibbl, 9)
      isto=lbbl(ibbl,10)
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
      if(lbbl(ibbl,1).ne.me) then
      if(lbbl(ibbl,11).gt.0) then
      jwfr=lbbl(ibbl,5)
      jwto=jwfr+lbbl(ibbl,12)-lbbl(ibbl,11)
      endif
      if(lbbl(ibbl,13).gt.0) then
      jsfr=lbbl(ibbl,6)
      jsto=jsfr+lbbl(ibbl,14)-lbbl(ibbl,13)
      endif
      endif
      lself=lbbl(ibbl,2).eq.lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
c
      handle_iw=lbbl(ibbl,15)
      handle_w=lbbl(ibbl,16)
      handle_is=lbbl(ibbl,17)
      handle_s=lbbl(ibbl,18)
      handle_pw=lbbl(ibbl,19)
      handle_ps=lbbl(ibbl,20)
      handle_pwp=lbbl(ibbl,21)
      handle_psp=lbbl(ibbl,22)
c
      if(handle_iw.ne.0.or.handle_w.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_iw)
      call ga_nbwait(handle_w)
      if(npackw.ne.0) then
      call argos_space_unpackw(jwto-jwfr+1,iwl(indexw,1),iwlp(indexw,1))
      endif
      endif
      if(handle_is.ne.0.or.handle_s.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_is)
      call ga_nbwait(handle_s)
      if(npackw.ne.0) then
      call argos_space_unpack(jsto-jsfr+1,isl(indexs,1),islp(indexs,1))
      endif
      endif
      if(handle_pw.ne.0.or.handle_ps.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_pw)
      call ga_nbwait(handle_ps)
      endif
      if(handle_pwp.ne.0.or.handle_psp.ne.0) then
      lnew=.true.
      call ga_nbwait(handle_pwp)
      call ga_nbwait(handle_psp)
      endif
c
      return
      end
c $Id$
