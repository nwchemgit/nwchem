      logical function grid_insidesp(nx,ny,nz,l_rect,rad)
c
c$Id: grid_sphere.F,v 1.9 2000-03-29 01:59:11 edo Exp $
c
      implicit none
      integer nx,ny,nz            ! [in]
      double precision l_rect,rad ! [in]
c
      double precision ndx,ndy,ndz,dista,nxyz,nxyz2
      grid_insidesp=.true.
      nxyz=nint(rad/l_rect)*2d0
      nxyz2=nxyz/2
      if(nx.le.nxyz2) then
         ndx=nxyz2-nx
      else
         ndx=dble(nx-nxyz2)-1d0
      endif
      if(ny.le.nxyz2) then
         ndy=nxyz2-ny
      else
         ndy=dble(ny-nxyz2)-1d0
      endif
      if(nz.le.nxyz2) then
         ndz=nxyz2-nz
      else
         ndz=dble(nz-nxyz2)-1d0
      endif
      dista=l_rect*sqrt(ndx*ndx+ndy*ndy+ndz*ndz)
      if(dista.gt.rad) grid_insidesp=.false.
      return
      end
      integer function grid_ncubes(rad,ynxyz,fnxyz,l_rect)
      implicit none
      integer nx,ny,nz,nxyz
      integer ninside
      double precision rad ! [in]
      integer fnxyz(*)     ! [out]
      logical ynxyz        ! [in] get points
      logical grid_insidesp
      double precision l_rect ! [in]
      external grid_insidesp
#include "cgridfile.fh"
#include "cdft.fh"
      integer nxyz2n,onx,ony,onz,onxyz
      nxyz2n(onx,ony,onz,onxyz)=onx-1+
     +(ony-1+(onz-1)*onxyz)*onxyz
      nxyz=nint(rad/l_rect)*2d0
      ninside=0
      do nx=1,nxyz
         do ny=1,nxyz
            do nz=1,nxyz
               if(grid_insidesp(nx,ny,nz,l_rect,rad)) then
                  ninside=ninside+1
                  if(ynxyz) then
                     fnxyz(ninside)=nxyz2n(nx,ny,nz,nxyz)
                     call grid_n2nxyz(nxyz,fnxyz(ninside),onx,ony,onz)
                     if(onx.ne.nx.or.ony.ne.ny.or.onz.ne.nz) then
                        write(0,*) ' IN nx ny nz ',nx,ny,nz,
     ,                    fnxyz(ninside)
                        write(0,*) ' OUT nx ny nz ',onx,ony,onz
                        call errquit('grid_ncube: mapping wrong',9)
                     endif
                  endif
              endif
            enddo
         enddo
      enddo
!      write(6,*) ' nxyz ninside',nxyz,ninside
      grid_ncubes=ninside
      return
      end
      subroutine grid_n2nxyz(nxyz,npoint,nx,ny,nz)
      implicit none
      integer nx,ny,nz    ! [out]
      integer nxyz,npoint ! [in]
c
      integer nresid
c      
      nz=npoint/(nxyz*nxyz)+1
      nresid=npoint-(nz-1)*nxyz*nxyz
      ny=nresid/nxyz+1
      nresid=nresid-(ny-1)*nxyz
      nx=nresid+1
      return
      end
      subroutine grid_rmaxmin(nx,ny,nz,rmax,rmin,l,nxyz)
      implicit none
      integer nx,ny,nz      ! [in]
      integer nxyz ! [in]
      double precision l    ! [in]
      double precision rmax ! [out]
      double precision rmin ! [out]
c
      double precision onx,ony,onz,nxyz2
      double precision onxyz2
c
      onx=nx
      ony=ny
      onz=nz
      rmax=0d0
      rmin=0d0
      nxyz2=nxyz/2
      onxyz2=nxyz*0.5d0
      if(nx.le.nxyz2) then
         rmax=rmax+((onxyz2-onx)+1d0)**2
         rmin=rmin+(onxyz2-onx)**2
      elseif(nx.gt.nxyz2) then
         rmax=rmax+(onx-onxyz2)**2
         rmin=rmin+(onx-onxyz2-1d0)**2
      endif
      if(ny.le.nxyz2) then
         rmax=rmax+((onxyz2-ony)+1d0)**2
         rmin=rmin+(onxyz2-ony)**2
      elseif(ny.gt.nxyz2) then
         rmax=rmax+(ony-onxyz2)**2
         rmin=rmin+(ony-onxyz2-1d0)**2
      endif
      if(nz.le.nxyz2) then
         rmax=rmax+((onxyz2-onz)+1d0)**2
         rmin=rmin+(onxyz2-onz)**2
      elseif(nz.gt.nxyz2) then
         rmax=rmax+(onz-onxyz2)**2
         rmin=rmin+(onz-onxyz2-1d0)**2
      endif

      rmax=sqrt(rmax)*l
      rmin=sqrt(rmin)*l
      return
      end
      logical function grid_incube(x,y,z,nx,ny,nz,nxyz,ll)
      double precision x,y,z ! [in] point coors
      integer nx,ny,nz       ! [in] cube definition
      double precision ll ! [in]
c
c     returns true if xyz is inside cube defined by nx,ny,nz
c
#include "cdft.fh"      
c     
      integer nxyz,nxyz2
c
      nxyz2=nxyz/2
      grid_incube=.false.
      if(x.le.(nx-nxyz2)*ll.and.x.gt.(nx-nxyz2-1)*ll) then
         if(y.le.(ny-nxyz2)*ll.and.y.gt.(ny-nxyz2-1)*ll) then
            if(z.le.(nz-nxyz2)*ll.and.z.gt.(nz-nxyz2-1)*ll) then
               grid_incube=.true.
            endif
         endif
      endif
c      write(0,1) x,y,z,nx,ny,nz,grid_incube
c 1    format(3f20.12,3I3,I5)
      return
      end
      double precision function grid_lcube(rad)
      double precision rad
c
      double precision l_cube
      parameter(l_cube=5.5d0)
      integer n
c
      n=(rad+rad)/l_cube
      grid_lcube=dble(rad+rad)/dble(n)
      return
      end
