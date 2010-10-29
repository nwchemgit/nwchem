      subroutine hnd_spdgia(ng,nr,g,xyz)
c
c $Id$
c
      implicit none
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
c
c     Input parameters
c
      integer ng,nr
      double precision g(ngint,6),xyz(ng*nr,3,3,*)
c
c     Local variables
c
      integer nn,nx,ny,nz,igr,Nxyz(3)
      double precision dum1,dum2,dum3,dum4,dum5,dum6,dum7,fij,fkl
      integer llkjt,llkt,ni,nj,nk,nl,i,j,k,l,m
      integer ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
c
c     ----- sum up ( ix * iy * iz ) -----
c
      llkt   = lmt*(lkt+1)
      llkjt  = llkt*ljt
      nk=numl             
      nj=nk*numk
      ni=nj*numj
      ijkl=0
      do 240 i=1,numi
        call getNxyz(lit-1,i,Nxyz)
        ix = Nxyz(1)
        iy = Nxyz(2)
        iz = Nxyz(3)
        do 230 j=1,numj
          call getNxyz(ljt-1,j,Nxyz)
          jx = Nxyz(1)
          jy = Nxyz(2)
          jz = Nxyz(3)
          do 220 k=1,numk
            call getNxyz(lkt-1,k,Nxyz)
            kx = Nxyz(1)
            ky = Nxyz(2)
            kz = Nxyz(3)
            do 210 l=1,numl
              call getNxyz(lmt-1,l,Nxyz)
              lx = Nxyz(1)
              ly = Nxyz(2)
              lz = Nxyz(3)
              ijkl=ijkl+1
              nn=1 + ni*(i-1) + nj*(j-1) + nk*(k-1) + (l-1)
              nx=1 + ix*llkjt + jx*llkt  + kx*lmt   + lx
              ny=1 + iy*llkjt + jy*llkt  + ky*lmt   + ly
              nz=1 + iz*llkjt + jz*llkt  + kz*lmt   + lz
              dum1=0.0d0
              dum2=0.0d0
              dum3=0.0d0
              dum4=0.0d0
              dum5=0.0d0
              dum6=0.0d0
              dum7=0.0d0
              do igr=1,ng*nr
               dum1=dum1+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
               dum2=dum2+xyz(igr,1,2,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
               dum3=dum3+xyz(igr,1,1,nx)*xyz(igr,2,2,ny)*xyz(igr,3,1,nz)
               dum4=dum4+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,2,nz)
               dum5=dum5+xyz(igr,1,3,nx)*xyz(igr,2,1,ny)*xyz(igr,3,1,nz)
               dum6=dum6+xyz(igr,1,1,nx)*xyz(igr,2,3,ny)*xyz(igr,3,1,nz)
               dum7=dum7+xyz(igr,1,1,nx)*xyz(igr,2,1,ny)*xyz(igr,3,3,nz)
              enddo     
              g(nn,1)=g(nn,1)+(qijx*dum1+tijy*dum4-tijz*dum3)
              g(nn,2)=g(nn,2)+(qijy*dum1+tijz*dum2-tijx*dum4)
              g(nn,3)=g(nn,3)+(qijz*dum1+tijx*dum3-tijy*dum2)
              g(nn,4)=g(nn,4)+(qklx*dum1+tkly*dum7-tklz*dum6)
              g(nn,5)=g(nn,5)+(qkly*dum1+tklz*dum5-tklx*dum7)
              g(nn,6)=g(nn,6)+(qklz*dum1+tklx*dum6-tkly*dum5)
  210       continue
  220     continue
  230   continue
  240 continue
      return
      end
