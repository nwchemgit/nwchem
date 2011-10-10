      subroutine argos_space_ownerw(xw,iwl,numwm,boxsiz,ibownr)
c
      implicit none
c
#include "argos_space_common.fh"
c
      real*8 xw(mwm,3,mwa)
      real*8 boxsiz(maxbox,3)
      integer iwl(mwm,miw2),ibownr(maxbox,3)
      integer numwm
c
      integer iwm,ibx,iby,ibz,i,k,ipx,ipy,ipz
      real*8 cgx,cgy,cgz,xt(3)
c
c     determine the box and owning node for each solvent molecule
c
      do 1 iwm=1,numwm
      cgx=zero
      cgy=zero
      cgz=zero
      do 2 k=1,nwa
      cgx=cgx+xw(iwm,1,k)
      cgy=cgy+xw(iwm,2,k)
      cgz=cgz+xw(iwm,3,k)
    2 continue
      ibx=0
      iby=0
      ibz=0
      if(nbxtyp.ne.1) then
      xt(1)=cgx
      xt(2)=cgy
      xt(3)=cgz
      else
      xt(1)=box(1)*(vlati(1,1)*cgx+vlati(1,2)*cgy+vlati(1,3)*cgz)
      xt(2)=box(2)*(vlati(2,1)*cgx+vlati(2,2)*cgy+vlati(2,3)*cgz)
      xt(3)=box(3)*(vlati(3,1)*cgx+vlati(3,2)*cgy+vlati(3,3)*cgz)
      endif
      do 3 i=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(i,1)) ibx=i
    3 continue
      do 4 i=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(i,2)) iby=i
    4 continue
      do 5 i=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(i,3)) ibz=i
    5 continue
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
      iwl(iwm,lwbox)=(ibz*nby+iby)*nbx+ibx
      iwl(iwm,lwnod)=(ipz*npy+ipy)*npx+ipx
    1 continue
c
      return
      end
c $Id$
