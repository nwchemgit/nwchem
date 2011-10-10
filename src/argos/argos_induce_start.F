      subroutine argos_induce_start(iwl,isl,xw,xwm,xs,pw,pwp,ps,psp)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
c
      logical lself,lpbcs
      integer ibbl,nbbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer i,j,k
      real*8 pmax
c
c     initialize induced fields to zero for first order polarization
c
      if(ipolt.eq.1) then
      if(nwmloc.gt.0) then
      do 1 k=1,nwa
      do 2 j=1,3
      do 3 i=1,nwmloc
      pw(i,j,k,1)=zero
      pw(i,j,k,2)=zero
    3 continue
    2 continue
    1 continue
      endif
      if(nsaloc.gt.0) then
      do 4 j=1,3
      do 5 i=1,nsaloc
      ps(i,j,1)=zero
      ps(i,j,2)=zero
    5 continue
    4 continue
      endif
      if(lpert2.or.lpert3) then
      if(nwmloc.gt.0) then
      do 6 k=1,nwa
      do 7 j=1,3
      do 8 i=1,nwmloc
      pwp(i,j,k,1,1)=zero
      pwp(i,j,k,1,2)=zero
      pwp(i,j,k,2,1)=zero
      pwp(i,j,k,2,2)=zero
    8 continue
    7 continue
    6 continue
      endif
      if(nsaloc.gt.0) then
      do 9 j=1,3
      do 10 i=1,nsaloc
      psp(i,j,1,1)=zero
      psp(i,j,1,2)=zero
      psp(i,j,2,1)=zero
      psp(i,j,2,2)=zero
   10 continue
    9 continue
      endif
      endif
      endif
c
c     iterative cycle to generate induced fields
c     ------------------------------------------
c
      npolit=0
      call argos_space_nbbl(nbbl)
   11 continue
      npolit=npolit+1
c
c     copy fields from previous iteration
c     -----------------------------------
c
      if(nwmloc.gt.0) then
      do 12 k=1,nwa
      do 13 j=1,3
      do 14 i=1,nwmloc
      pw(i,j,k,2)=pw(i,j,k,1)
      pw(i,j,k,1)=zero
   14 continue
   13 continue
   12 continue
      endif
      if(nsaloc.gt.0) then
      do 15 j=1,3
      do 16 i=1,nsaloc
      ps(i,j,2)=ps(i,j,1)
      ps(i,j,1)=zero
   16 continue
   15 continue
      endif
      if(mdtype.gt.3) then
      if(nwmloc.gt.0) then
      do 17 k=1,nwa
      do 18 j=1,3
      do 19 i=1,nwmloc
      pwp(i,j,k,1,2)=pwp(i,j,k,1,1)
      pwp(i,j,k,1,1)=zero
      pwp(i,j,k,2,2)=pwp(i,j,k,2,1)
      pwp(i,j,k,2,1)=zero
   19 continue
   18 continue
   17 continue
      endif
      if(nsaloc.gt.0) then
      do 20 j=1,3
      do 21 i=1,nsaloc
      psp(i,j,1,2)=psp(i,j,1,1)
      psp(i,j,1,1)=zero
      psp(i,j,2,2)=psp(i,j,2,1)
      psp(i,j,2,1)=zero
   21 continue
   20 continue
      endif
      endif
c
c     copy current fields into local global array
c     -------------------------------------------
c
      call argos_space_putp(me,pw,pwp,nwmloc,ps,psp,nsaloc,
     + lpert2.or.lpert3)
c
c     synchronize to ensure induced fields are available
c     --------------------------------------------------
c
      call ga_sync()
c
      do 22 ibbl=1,nbbl
c
      call argos_space_getxpbl(ibbl,lhandl,
     + iwl,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,xs,ps,psp,isfr,isto,jsfr,jsto,lself,lpbcs)
c
      call argos_cafe_comw(xw,xwm,jwfr,jwto)
c
      call argos_induce(lself,lpbcs,xw,xwm,pw,pwp,iwl(1,lwdyn),
     + int_mb(i_iwz),iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),ps,psp,
     + isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),isl(1,lsfrc),
     + isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),isl(1,lsct3),
     + int_mb(i_isz),isfr,isto,jsfr,jsto,lpbc,lhandl)
c
      call argos_space_accpbl(ibbl,lhandl,pw,pwp,ps,psp,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
c
   22 continue
      lpair=.false.
      lload=.false.
c
      if(np.gt.0) call ga_sync()
c
      if(nwmloc.gt.0) then
      do 26 k=1,nwa
      do 27 j=1,3
      do 28 i=1,nwmloc
      pw(i,j,k,2)=pw(i,j,k,1)
      pw(i,j,k,1)=zero
   28 continue
   27 continue
   26 continue
      endif
      if(nsaloc.gt.0) then
      do 29 j=1,3
      do 30 i=1,nsaloc
      ps(i,j,2)=ps(i,j,1)
      ps(i,j,1)=zero
   30 continue
   29 continue
      endif
c
      call argos_space_getp(me,pw,pwp,nwmloc,ps,psp,nsaloc,
     + lpert2.or.lpert3,1)
c
      if(nwmloc.gt.0) then
      do 31 k=1,nwa
      do 32 j=1,3
      do 33 i=1,nwmloc
      pw(i,j,k,1)=pw(i,j,k,1)+pw(i,j,k,2)
      pw(i,j,k,2)=zero
   33 continue
   32 continue
   31 continue
      endif
      if(nsaloc.gt.0) then
      do 34 j=1,3
      do 35 i=1,nsaloc
      ps(i,j,1)=ps(i,j,1)+ps(i,j,2)
      ps(i,j,2)=zero
   35 continue
   34 continue
      endif
c
      call argos_space_getp(me,pw,pwp,nwmloc,ps,psp,nsaloc,
     + lpert2.or.lpert3,2)
c
      pmax=0.0d0
      do 23 k=1,nwa
      do 24 j=1,3
      do 25 i=1,nwmloc
      pmax=max(pmax,abs(pw(i,j,k,2)-pw(i,j,k,1)))
   25 continue
   24 continue
   23 continue
c
      if(np.gt.1) call ga_dgop(mrg_d06,pmax,1,'max')
c
      if(pmax.gt.ptol.and.npolit.le.mpolit.and.ipolt.gt.1) goto 11
c
c     copy current fields into local global array
c     -------------------------------------------
c
      call argos_space_putp(me,pw,pwp,nwmloc,ps,psp,nsaloc,
     + lpert2.or.lpert3)
c
      return
      end
c $Id$
