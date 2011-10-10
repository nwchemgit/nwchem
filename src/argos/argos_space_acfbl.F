      subroutine argos_space_acfbl(ibbl,lhandl,fw,fs,lbbl,jpl,
     + lpair,iwz,isz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer lbbl(mbbl,mbb2)
      integer jpl(mbox,mip2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,j,k,l,m,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      integer jwfr,jwto,jsfr,jsto
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call md_abort('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lbbl(ibbl,4)=lhandl
c
      if(jnode.eq.me) return
      if(ibbl.lt.nbbl) then
      if(lbbl(ibbl+1,1).eq.jnode.and.lbbl(ibbl+1,2).eq.jbox) return
      endif
c
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0)
     + call md_abort('Remote processor cell not found (2) ',0)
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+6*mwa+3,jl+9*mwa+2,
     + fw(nwmloc+1,1,1,1),mwm,one)
      if(llong) call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
     + jl+12*mwa+2,fw(nwmloc+1,1,1,2),mwm,one)
      if(lpair) then
      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
      call ga_acc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,iwz(nwmloc+1),mwm,1)
      endif
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (5)',nsaloc+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+6,jl+8,
     + fs(nsaloc+1,1,1),msa,one)
      if(llong) call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+9,jl+11,
     + fs(nsaloc+1,1,2),msa,one)
      if(lpair) then
      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
      call ga_acc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,isz(nsaloc+1),msa,1)
      endif
      endif
c
      m=1
      if(llong) m=2
c
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=nwmloc+1,mwm
      fw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
c
      do 7 j=1,3
      do 8 i=nsaloc+1,msa
      fs(i,j,l)=zero
    8 continue
    7 continue
    9 continue
c
      do 10 i=nwmloc+1,mwm
      iwz(i)=0
   10 continue
      do 11 i=nsaloc+1,msa
      isz(i)=0
   11 continue
c
      return
      end
c $Id$
