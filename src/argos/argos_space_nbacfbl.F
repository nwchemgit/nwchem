      subroutine argos_space_nbacfbl(ibbl,lhandl,fw,fs,lbbl,jpl,
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
      integer jwfr,jwto,jsfr,jsto,indexw,indexs
      integer handle_iw,handle_fw,handle_is,handle_fs
      integer handle_fwl,handle_fsl
c
      handle_iw=0
      handle_is=0
      handle_fw=0
      handle_fs=0
      handle_fwl=0
      handle_fsl=0
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
c      ndx=0
c      do 3 i=1,jpl(1,1)
c      if(jpl(30+i,1).eq.jbox) ndx=30+i
c    3 continue
c      if(ndx.eq.0)
c     + call md_abort('Remote processor cell not found (3) ',0)
c
      jwfr=lbbl(ibbl,11)
      jwto=lbbl(ibbl,12)
      jsfr=lbbl(ibbl,13)
      jsto=lbbl(ibbl,14)
c
c      if(lbbl(ibbl,1).ne.me) then
c      if(lbbl(ibbl,12).gt.0) then
c      jwfr=lbbl(ibbl,5)
c      jwto=jwfr+lbbl(ibbl,12)-lbbl(ibbl,11)
c      endif
c      if(lbbl(ibbl,14).gt.0) then
c      jsfr=lbbl(ibbl,6)
c      jsto=jsfr+lbbl(ibbl,14)-lbbl(ibbl,13)
c      endif
c      endif
c
      indexw=lbbl(ibbl,5)
      indexs=lbbl(ibbl,6)
c
      m=1
      if(llong) m=2
c
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(indexw+nwnon.gt.mwm)
     + call md_abort('Dimension mwm too small',indexw+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_w,il+jwfr-1,il+jwto-1,jl+6*mwa+3,jl+9*mwa+2,
     + fw(indexw,1,1,1),mwm,one,handle_fw)
      if(llong) call ga_nbacc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
     + jl+12*mwa+2,fw(indexw,1,1,2),mwm,one,handle_fwl)
      if(lpair) then
      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,
     + iwz(indexw),mwm,1,handle_iw)
      endif
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=indexw,indexw+nwnon-1
      fw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
    9 continue
      do 10 i=indexw,indexw+nwnon-1
      iwz(i)=0
   10 continue
      endif
c
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(indexs+nsnon.gt.msa)
     + call md_abort('Dimension msa too small (5)',indexs+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_s,il+jsfr-1,il+jsto-1,jl+6,jl+8,
     + fs(indexs,1,1),msa,one,handle_fs)
      if(llong) call ga_nbacc(ga_s,il+jsfr-1,il+jsto-1,jl+9,jl+11,
     + fs(indexs,1,2),msa,one,handle_fsl)
      if(lpair) then
      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
      call ga_nbacc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,
     + isz(indexs),msa,1,handle_is)
      endif
      do 12 l=1,m
      do 7 j=1,3
      do 8 i=indexs,indexs+nsnon-1
      fs(i,j,l)=zero
    8 continue
    7 continue
   12 continue
      do 11 i=indexs,indexs+nsnon-1
      isz(i)=0
   11 continue
      endif
c
      lbbl(ibbl,15)=handle_iw
      lbbl(ibbl,16)=handle_fw
      lbbl(ibbl,17)=handle_is
      lbbl(ibbl,18)=handle_fs
      lbbl(ibbl,19)=handle_fwl
      lbbl(ibbl,20)=handle_fsl
c
      return
      end
c $Id$
