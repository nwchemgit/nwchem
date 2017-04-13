      subroutine argos_space_owners(xs,isl,numsa,boxsiz,ibownr)
c
      implicit none
c
#include "argos_space_common.fh"
c
c
      real*8 xs(msa,3)
      real*8 boxsiz(maxbox,3)
      integer ibownr(maxbox,3)
      integer isl(msa,mis2)
      integer numsa,nfold
c
      integer isa,isaf,isal,jsa,ibox,inod,k
      integer ibx,iby,ibz,i,ipx,ipy,ipz
      real*8 xcgx,xcgy,xcgz,boxi(3),xcg(3)
      real*8 xcgmax,xcgmin,factor
c
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
      nfold=0
c
      isaf=1
      isal=0
      do 1 isa=1,numsa
c
      if(isa.lt.numsa) then
      if(isl(isa+1,lssgm).ne.isl(isaf,lssgm)) isal=isa
      else
      isal=isa
      endif
c
      if(isal.gt.0) then
c
      do 322 k=1,3
      xcgmin=xs(isaf,k)
      xcgmax=xs(isaf,k)
      do 323 jsa=isaf,isal
      xcgmin=min(xcgmin,xs(jsa,k))
      xcgmax=max(xcgmax,xs(jsa,k))
  323 continue
      xcg(k)=0.5d0*(xcgmax+xcgmin)
  322 continue
c
c      xcgx=xcg(1)
c      xcgy=xcg(2)
c      xcgz=xcg(3)
      xcgx=zero
      xcgy=zero
      xcgz=zero
      do 2 jsa=isaf,isal
      xcgx=xcgx+xs(jsa,1)
      xcgy=xcgy+xs(jsa,2)
      xcgz=xcgz+xs(jsa,3)
    2 continue
      factor=one/dble(isal-isaf+1)
      xcgx=factor*xcgx
      xcgy=factor*xcgy
      xcgz=factor*xcgz
cx
      xcgx=xcg(1)
      xcgy=xcg(2)
      xcgz=xcg(3)
c
      if(npbtyp.ne.0) then
      if(abs(xcgx).gt.boxh(1)) then
      xcgx=xcgx-nint(xcgx*boxi(1))*box(1)
      nfold=1
      endif
      if(abs(xcgy).gt.boxh(2)) then
      xcgy=xcgy-nint(xcgy*boxi(2))*box(2)
      nfold=1
      endif
      if(abs(xcgz).gt.boxh(3)) then
      xcgz=xcgz-nint(xcgz*boxi(3))*box(3)
      nfold=1
      endif
      endif
c
      ibx=0
      iby=0
      ibz=0
      do 3 i=1,nbx-1
      if(xcgx+boxh(1).gt.boxsiz(i,1)) ibx=i
    3 continue
      do 4 i=1,nby-1
      if(xcgy+boxh(2).gt.boxsiz(i,2)) iby=i
    4 continue
      do 5 i=1,nbz-1
      if(xcgz+boxh(3).gt.boxsiz(i,3)) ibz=i
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
      ibox=(ibz*nby+iby)*nbx+ibx
      inod=(ipz*npy+ipy)*npx+ipx
      do 6 jsa=isaf,isal
      isl(jsa,lsbox)=ibox
      isl(jsa,lsnod)=inod
    6 continue
c
      isaf=isa+1
      isal=0
c
      endif
c
    1 continue
c
      return
      end
c $Id$
