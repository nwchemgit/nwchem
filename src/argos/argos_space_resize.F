      subroutine argos_space_resize(ipmin,ibownr,boxsiz,
     + factor,facx,facy,facz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ipmin,ibownr(maxbox,3)
      real*8 factor,facx,facy,facz,boxsiz(maxbox,3)
c
      integer ibx,iby,ibz,ipminx,ipminy,ipminz
      real*8 bxi,bxr,byi,byr,bzi,bzr,rfact
c
c     load balancing based on resizing of node domains
c     ------------------------------------------------
c
      ipminx=mod(ipmin,npx)
      ipminy=mod((ipmin-ipminx)/npx,npy)
      ipminz=((ipmin-ipminx)/npx-ipminy)/npy
c
c     resizing boxes in x-dimension
c
      if(npx.gt.1.and.isload.eq.0) then
      do 7 ibx=nbx,2,-1
      boxsiz(ibx,1)=boxsiz(ibx,1)-boxsiz(ibx-1,1)
    7 continue
      bxi=zero
      bxr=zero
      do 8 ibx=1,nbx
      if(ibownr(ibx,1).eq.ipminx) then
      facx=1.0d0/boxsiz(ibx,1)
      boxsiz(ibx,1)=factor*boxsiz(ibx,1)
      if(boxsiz(ibx,1).lt.bxmin) boxsiz(ibx,1)=bxmin
      facx=facx*boxsiz(ibx,1)
      bxi=bxi+boxsiz(ibx,1)
      else
      bxr=bxr+boxsiz(ibx,1)
      endif
    8 continue
      rfact=(box(1)-bxi)/bxr
      do 9 ibx=1,nbx
      if(ibownr(ibx,1).ne.ipminx) boxsiz(ibx,1)=rfact*boxsiz(ibx,1)
      if(ibx.gt.1) boxsiz(ibx,1)=boxsiz(ibx,1)+boxsiz(ibx-1,1)
    9 continue
      endif
c
c     resizing boxes in y-dimension
c
      if(npy.gt.1.and.isload.eq.0) then
      do 10 iby=nby,2,-1
      boxsiz(iby,2)=boxsiz(iby,2)-boxsiz(iby-1,2)
   10 continue
      byi=zero
      byr=zero
      do 11 iby=1,nby
      if(ibownr(iby,2).eq.ipminy) then
      facy=1.0d0/boxsiz(iby,2)
      boxsiz(iby,2)=factor*boxsiz(iby,2)
      if(boxsiz(iby,2).lt.bymin) boxsiz(iby,2)=bymin
      facy=facy*boxsiz(iby,2)
      byi=byi+boxsiz(iby,2)
      else
      byr=byr+boxsiz(iby,2)
      endif
   11 continue
      rfact=(box(2)-byi)/byr
      do 12 iby=1,nby
      if(ibownr(iby,2).ne.ipminy) boxsiz(iby,2)=rfact*boxsiz(iby,2)
      if(iby.gt.1) boxsiz(iby,2)=boxsiz(iby,2)+boxsiz(iby-1,2)
   12 continue
      endif
c
c     resizing boxes in z-dimension
c
      if(npz.gt.1) then
      do 13 ibz=nbz,2,-1
      boxsiz(ibz,3)=boxsiz(ibz,3)-boxsiz(ibz-1,3)
   13 continue
      bzi=zero
      bzr=zero
      do 14 ibz=1,nbz
      if(ibownr(ibz,3).eq.ipminz) then
      facz=1.0d0/boxsiz(ibz,3)
      boxsiz(ibz,3)=factor*boxsiz(ibz,3)
      if(boxsiz(ibz,3).lt.bzmin) boxsiz(ibz,3)=bzmin
      facz=facz*boxsiz(ibz,3)
      bzi=bzi+boxsiz(ibz,3)
      else
      bzr=bzr+boxsiz(ibz,3)
      endif
   14 continue
      rfact=(box(3)-bzi)/bzr
      do 15 ibz=1,nbz
      if(ibownr(ibz,3).ne.ipminz) boxsiz(ibz,3)=rfact*boxsiz(ibz,3)
      if(ibz.gt.1) boxsiz(ibz,3)=boxsiz(ibz,3)+boxsiz(ibz-1,3)
   15 continue
      endif
c
      call ga_brdcst
     + (msp_19,boxsiz,ma_sizeof(mt_dbl,3*maxbox,mt_byte),0)
c
      return
      end
c $Id$
