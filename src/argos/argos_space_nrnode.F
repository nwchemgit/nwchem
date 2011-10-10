      subroutine argos_space_nrnode
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer ix,iy,iz,npt,i,j,k,l
c
c     this routine distributes the available processes in the
c     Cartesian directions
c
c     npx : number of processes in x direction
c     npy : number of processes in y direction
c     npz : number of processes in z direction
c
c     determine node dimensions
c
      if(npx*npy*npz.ne.np) then
      if(npx+npy+npz.gt.0)
     + call md_abort('Specified npx*npy*npz ne np',0)
c
      npt=0
      do 1 i=1,np
      do 2 j=i,np
      do 3 k=j,np
      l=i*j*k
      if(l.eq.np) then
      if(l.eq.npt) then
      if(k.gt.npz) goto 3
      if(i+j+k.lt.npx+npy+npz) then
      npx=i
      npy=j
      npz=k
      endif
      goto 3
      else
      npt=np
      npx=i
      npy=j
      npz=k
      endif
      endif
    3 continue
    2 continue
    1 continue
c
      if(npx*npy*npz.ne.np) call md_abort('nrnode: code error',0)
      endif
c
c     determine processor location of me
c
      do 4 ix=1,npx
      do 5 iy=1,npy
      do 6 iz=1,npz
      if(me.eq.((iz-1)*npy+(iy-1))*npx+(ix-1)) then
      mex=ix-1
      mey=iy-1
      mez=iz-1
      endif
    6 continue
    5 continue
    4 continue
c
      return
      end
c $Id$
