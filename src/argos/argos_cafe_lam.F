      subroutine argos_cafe_lam(lamtyp,lambda,maxlam,explam,lfnout,
     + rlambd,dlambd,
     + wgt,vdw,chg,mwb,nwb,nbp,bndw,
     + mwh,nwh,nhp,angw,mwd,nwd,ndp,dihw,mwo,nwo,nop,oopw,
     + msb,npb,bnds,msh,nph,angs,msd,npd,dihs,mso,npo,oops,
     + lfnpmf)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lamtyp,lambda,maxlam,lfnout
      real*8 explam,rlambd,dlambd
      real*8 wgt(mat,mset),vdw(mat,mat,map,mset),chg(mqt,mqp,mset)
      integer mwb,mwh,mwd,mwo,nbp,nhp,ndp,nop
      integer nwb,nwh,nwd,nwo
      integer msb,msh,msd,mso,npb,nph,npd,npo
      real*8 bndw(mwb,nbp,6),angw(mwh,nhp,6)
      real*8 dihw(mwd,ndp,6),oopw(mwo,nop,6)
      real*8 bnds(msb,npb,6),angs(msh,nph,6)
      real*8 dihs(msd,npd,6),oops(mso,npo,6)
      integer lfnpmf
c
      real*8 tlam(3),rlam(3),slam(3)
      integer i,j,k,l
      character*10 pdate,ptime
c
      tlam(1)=dble(lambda-1)/dble(maxlam-1)
      tlam(2)=(dble(lambda-1)-half)/dble(maxlam-1)
      tlam(3)=(dble(lambda-1)+half)/dble(maxlam-1)
      if(tlam(2).lt.zero) tlam(2)=zero
      if(tlam(3).gt.one) tlam(3)=one
c
      if(iabs(lamtyp).eq.1) then
      rlam(1)=tlam(1)
      rlam(2)=tlam(2)
      rlam(3)=tlam(3)
      elseif(iabs(lamtyp).eq.2) then
      rlam(1)=tlam(1)**explam
      rlam(2)=tlam(2)**explam
      rlam(3)=tlam(3)**explam
      elseif(iabs(lamtyp).eq.3) then
      rlam(1)=one-(one-tlam(1))**explam
      rlam(2)=one-(one-tlam(2))**explam
      rlam(3)=one-(one-tlam(3))**explam
      endif
      if(lamtyp.lt.0) then
      rlam(1)=one-rlam(1)
      rlam(2)=one-rlam(2)
      rlam(3)=one-rlam(3)
      endif
      slam(1)=one-rlam(1)
      slam(2)=one-rlam(2)
      slam(3)=one-rlam(3)
c
      do 1 k=1,3
      do 2 i=1,mat
      wgt(i,k)=rlam(k)*wgt(i,6)+slam(k)*wgt(i,5)
      do 3 l=1,map
      do 4 j=1,mat
      vdw(i,j,l,k)=rlam(k)*vdw(i,j,l,6)+slam(k)*vdw(i,j,l,5)
    4 continue
    3 continue
    2 continue
      do 5 j=1,mqp
      do 6 i=1,mqt
      chg(i,j,k)=rlam(k)*chg(i,j,6)+slam(k)*chg(i,j,5)
    6 continue
    5 continue
      do 7 j=1,nbp
      do 8 i=1,mwb
      bndw(i,j,k)=rlam(k)*bndw(i,j,6)+slam(k)*bndw(i,j,5)
    8 continue
    7 continue
      do 9 j=1,nhp
      do 10 i=1,mwh
      angw(i,j,k)=rlam(k)*angw(i,j,6)+slam(k)*angw(i,j,5)
   10 continue
    9 continue
      do 11 j=1,ndp
      do 12 i=1,mwd
      dihw(i,j,k)=rlam(k)*dihw(i,j,6)+slam(k)*dihw(i,j,5)
   12 continue
   11 continue
      do 13 j=1,nop
      do 14 i=1,mwo
      oopw(i,j,k)=rlam(k)*oopw(i,j,6)+slam(k)*oopw(i,j,5)
   14 continue
   13 continue
      do 15 j=1,nbp
      do 16 i=1,msb
      bnds(i,j,k)=rlam(k)*bnds(i,j,6)+slam(k)*bnds(i,j,5)
   16 continue
   15 continue
      do 17 j=1,nhp
      do 18 i=1,msh
      angs(i,j,k)=rlam(k)*angs(i,j,6)+slam(k)*angs(i,j,5)
   18 continue
   17 continue
      do 19 j=1,ndp
      do 20 i=1,msd
      dihs(i,j,k)=rlam(k)*dihs(i,j,6)+slam(k)*dihs(i,j,5)
   20 continue
   19 continue
      do 21 j=1,nop
      do 22 i=1,mso
      oops(i,j,k)=rlam(k)*oops(i,j,6)+slam(k)*oops(i,j,5)
   22 continue
   21 continue
      shift0(k)=rlam(k)*shift0(6)+slam(k)*shift0(5)
      shift1(k)=rlam(k)*shift1(6)+slam(k)*shift1(5)
    1 continue
c
      if(lpmf) call argos_cafe_lampmf(rlam,slam,dbl_mb(i_rpmf),lfnpmf)
c
      if(me.eq.0) then
      call swatch(pdate,ptime)
      write(lfnout,1000) rlam(1),ith,pdate,ptime
 1000 format(/' THERMODYNAMIC INTEGRATION, LAMBDA=',f10.5,
     + 5x,12l1,1x,12l1,t110,2a10)
      endif
c
      rlambd=rlam(1)
      dlambd=rlam(3)-rlam(2)
c
      return
      end
c $Id$
