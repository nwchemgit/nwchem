      subroutine argos_space_decomp(ibownr,boxsiz,iburen,ibindx)
c
      implicit none
c
#include "argos_space_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
      integer iburen(np,27,2)
      integer ibindx(np)
c
      integer ibx,iby,ibz,i,j,ix,jx,iy,jy,iz,jz,jnode,nrnod,nbiown
c
c     check dimensions of ibownr
c
      if(nbx.gt.maxbox.or.nby.gt.maxbox.or.nbz.gt.maxbox)
     + call md_abort('Dimension maxbox too small',0)
c
c     determine the node dimension for each sub box
c
      do 1 ibx=1,nbx
      ibownr(ibx,1)=((ibx-1)*npx)/nbx
      boxsiz(ibx,1)=(box(1)*dble(ibx))/dble(nbx)
    1 continue
      do 2 iby=1,nby
      ibownr(iby,2)=((iby-1)*npy)/nby
      boxsiz(iby,2)=(box(2)*dble(iby))/dble(nby)
    2 continue
      do 3 ibz=1,nbz
      ibownr(ibz,3)=((ibz-1)*npz)/nbz
      boxsiz(ibz,3)=(box(3)*dble(ibz))/dble(nbz)
    3 continue
      if(iand(idebug,1).eq.1) then
      write(lfndbg,8003) nbx,nby,nbz,maxbox,npx,npy,npz
 8003 format('ibownr in argos_space_decomp',7i5)
      write(lfndbg,8004) (ibownr(ibx,1),ibx=1,nbx)
 8004 format('ibownr x',/,(20i5))
      write(lfndbg,8005) (ibownr(iby,2),iby=1,nby)
 8005 format('ibownr y',/,(20i5))
      write(lfndbg,8006) (ibownr(ibz,3),ibz=1,nbz)
 8006 format('ibownr z',/,(20i5))
      write(lfndbg,8002)
 8002 format('boxlist')
      endif
c
c     determine neighboring nodes and store in neighb(27,2)
c     such that:
c
c     neighb(n,1) is the n-th neighbor of node me
c     neighb(n,2) is the node of which node me is the n-th neighbor
c
c     a value of -1 indicates that such node does not exist
c
      do 4 i=1,27
      neighb(i,1)=-1
      neighb(i,2)=-1
    4 continue
c
      do 5 ix=1,3
      jx=mex+ix-2
      if(npbtyp.gt.0) then
      if(npx.gt.2.and.jx.lt.0) jx=jx+npx
      if(npx.gt.2.and.jx.ge.npx) jx=jx-npx
      endif
      if(jx.ge.0.and.jx.lt.npx) then
      do 6 iy=1,3
      jy=mey+iy-2
      if(npbtyp.gt.0) then
      if(npy.gt.2.and.jy.lt.0) jy=jy+npy
      if(npy.gt.2.and.jy.ge.npy) jy=jy-npy
      endif
      if(jy.ge.0.and.jy.lt.npy) then
      do 7 iz=1,3
      jz=mez+iz-2
      if(npbtyp.eq.1) then
      if(npz.gt.2.and.jz.lt.0) jz=jz+npz
      if(npz.gt.2.and.jz.ge.npz) jz=jz-npz
      endif
      if(jz.ge.0.and.jz.lt.npz) then
      jnode=((jz*npy)+jy)*npx+jx
      neighb(3*(3*(ix-1)+(iy-1))+iz,1)=jnode
      neighb(3*(3*(3-ix)+(3-iy))+4-iz,2)=jnode
      endif
    7 continue
      endif
    6 continue
      endif
    5 continue
c
      nbiown=30
      do 8 ibx=1,nbx
      do 9 iby=1,nby
      do 10 ibz=1,nbz
      nrnod=(ibownr(ibz,3)*npy+ibownr(iby,2))*npx+ibownr(ibx,1)
      if(me.eq.nrnod) nbiown=nbiown+1
   10 continue
    9 continue
    8 continue
c
      do 11 j=1,27
      do 12 i=1,np
      iburen(i,j,1)=0
      iburen(i,j,2)=0
   12 continue
   11 continue
      do 13 i=1,np
      ibindx(i)=0
   13 continue
      do 14 j=1,27
      iburen(me+1,j,1)=neighb(j,1)
      iburen(me+1,j,2)=neighb(j,2)
      if(neighb(j,1).ge.0) ibindx(neighb(j,1)+1)=j
   14 continue
c
      call ga_igop(msp_27,iburen,np*54,'+')
c
c     broadcast the maximum number of sub-boxes per node to all nodes
c
      if(nbiown.gt.mbox) call md_abort('Error in mbox',0)
c
      return
      end
c $Id$
