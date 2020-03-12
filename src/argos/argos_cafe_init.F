      subroutine argos_cafe_init(stimei,lp,llng,bx,vl,vli,zwi,zs,
     + eww,esw,ess,fss,esa)
c     + vw,vwt,vs,vst,zwi,zs,eww,esw,ess)
c
c     in r*8 : bx(3) : box dimensions
c     in log : lp    : flag to force pairlist recalculation
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
c      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),vs(msa,3),vst(msa,3)
      real*8 bx(3),vl(3,3),vli(3,3),zwi(3,3,2),zs(msf,3,3,2)
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2)
      real*8 fss(msf,msf,3,2)
      real*8 esa(nsa,2),stimei
      logical lp,llng
c
      integer i,j,k,l,m
c
c     initialize pme flags
c
      stime=stimei
      if(nhops.gt.0) call argos_cafe_hopinit(int_mb(i_lsthop),
     + dbl_mb(i_timhop))
c
      if(ipme.gt.0) then
      call argos_pme_flag(0,0,me)
      call argos_pme_flag(1,0,me)
      endif
c
      lpair=lp
      llong=llng
c
      ntest=0
      m=1
      if(llong) m=2
c
      do 1 i=1,3
      box(i)=bx(i)
      boxh(i)=half*bx(i)
      do 2 j=1,3
      vlat(i,j)=vl(i,j)
      vlati(i,j)=vli(1,3)
    2 continue
    1 continue
      volume=box(1)*box(2)*box(3)
c
      do 8 l=1,m
      do 3 i=1,mpe
      eww(i,l)=zero
      do 4 j=1,msf
      esw(j,i,l)=zero
      do 5 k=1,msf
      ess(k,j,i,l)=zero
    5 continue
    4 continue
    3 continue
      do 16 k=1,3
      do 6 i=1,3
      zw(k,i,l)=zero
      zwi(k,i,l)=zero
      do 7 j=1,msf
      zs(j,k,i,l)=zero
    7 continue
    6 continue
   16 continue
      if(m.eq.2) then
      vpme(1)=zero
      vpme(2)=zero
      vpme(3)=zero
      vpme(4)=zero
      vpme(5)=zero
      vpme(6)=zero
      endif
      if(ntype.eq.3) then
      do 9 i=1,24
      deriv(i,l)=zero
    9 continue
      endif
      ep2(l)=zero
      ep3(l)=zero
      do 21 i=1,msf
      do 22 j=1,msf
      fss(i,j,1,l)=zero
      fss(i,j,2,l)=zero
      fss(i,j,3,l)=zero
   22 continue
   21 continue
    8 continue
      ep2(3)=zero
      ep3(3)=zero
c
      if(npener.gt.0) then
      do 10 i=1,nsa
      esa(i,1)=zero
      esa(i,2)=zero
   10 continue
      endif
c
      if(me.eq.0) eww(6,1)=dble(nwm)*ewc(iset)
c
      lpww=1
      lpsw=1
      lpss=1
      if(llong) then
      lpww=2
      lpsw=2
      lpss=2
      endif
c
      if(lpair) then
      llww=0
      lsww=0
      llsw=0
      lssw=0
      llss=0
      lsss=0
      nlda=0
      endif
c
      do 23 i=1,4*nhop
      dbl_mb(i_uda-1+i)=zero
   23 continue
c
      if(.not.llist) lpair=.true.
c
c     allocate memory for the scratch arrays
c
      if(.not.lscr) call argos_cafe_salloc()
c
c     allocate memory for the pairlists
c
      if(.not.llst) call argos_cafe_lalloc()
c
      ndxp=0
c
      call argos_cafe_binit(dbl_mb(i_rbnd(2)),numb(2),2)
      call argos_cafe_binit(dbl_mb(i_rang(2)),numh(2),2)
      call argos_cafe_binit(dbl_mb(i_rub(2)),numh(2),2)
      call argos_cafe_binit(dbl_mb(i_rdih(2)),numd(2),2)
      call argos_cafe_binit(dbl_mb(i_rimp(2)),numi(2),2)
c
      return
      end
c $Id$
