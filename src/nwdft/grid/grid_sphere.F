      logical function grid_insidesp(nx,ny,nz,l_rect,rad)
c
c$Id: grid_sphere.F,v 1.5 2000-03-18 07:14:22 edo Exp $
c
      implicit none
      integer nx,ny,nz            ! [in]
      double precision l_rect,rad ! [in]
c
      double precision ndx,ndy,ndz,dista
      grid_insidesp=.true.
      ndx=0d0
      if(nx.gt.0) then
         ndx=dble(nx)-0.5d0
      elseif(nx.lt.0) then
         ndx=dble(nx)+0.5d0
      endif
      ndy=0d0
      if(ny.gt.0) then
         ndy=dble(ny)-0.5d0
      elseif(ny.lt.0) then
         ndy=dble(ny)+0.5d0
      endif
      ndz=0d0
      if(nz.gt.0) then
         ndz=dble(nz)-0.5d0
      elseif(nz.lt.0) then
         ndz=dble(nz)+0.5d0
      endif
      dista=l_rect*sqrt(ndx*ndx+ndy*ndy+ndz*ndz)
      if(dista.gt.rad) grid_insidesp=.false.
      return
      end
      integer function grid_ncubes(rad,ynxyz,fnxyz,l_rect)
      implicit none
      integer nx,ny,nz,nxyz
      integer ninside
      double precision l_rect ! [in]
      double precision rad ! [in]
      integer fnxyz(*)     ! [out]
      logical ynxyz        ! [in] get points
      logical grid_insidesp
      external grid_insidesp
#include "cgridfile.fh"
#include "cdft.fh"
      integer nxyz2n,onx,ony,onz,onxyz
      nxyz2n(onx,ony,onz,onxyz)=(onx+onxyz)+
     +( (ony+onxyz)+(onz+onxyz)*(2*onxyz+1) )*(2*onxyz+1)
      nxyz=nint((rad*2d0-l_rect)/(l_rect*2d0))
      write(0,*) ' rad l_rect',rad,l_rect,nxyz
      ninside=0
      do nx=-nxyz,nxyz
         do ny=-nxyz,nxyz
            do nz=-nxyz,nxyz
               if(grid_insidesp(nx,ny,nz,l_rect,rad)) then
                  ninside=ninside+1
                  if(ynxyz) then
                     fnxyz(ninside)=nxyz2n(nx,ny,nz,nxyz)
                     call grid_n2nxyz(
     .                    nxyz,fnxyz(ninside),onx,ony,onz)
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
      write(6,*) ' nxyz ninside',nxyz,ninside
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
      nz=npoint/((2*nxyz+1)*(2*nxyz+1))-nxyz
      nresid=npoint-(nz+nxyz)*(2*nxyz+1)*(2*nxyz+1)
      ny=nresid/(2*nxyz+1)-nxyz
      nresid=nresid-(ny+nxyz)*(2*nxyz+1)
      nx=nresid-nxyz
      return
      end
      subroutine grid_rmaxmin(nx,ny,nz,rmax,rmin,l)
      implicit none
      integer nx,ny,nz      ! [in]
      double precision l    ! [in]
      double precision rmax ! [out]
      double precision rmin ! [out]
c
      double precision rr,onx,ony,onz
c
      onx=nx
      ony=ny
      onz=nz
      
      rmax=sqrt( (abs(onx)+.5d0)**2d0+
     +(abs(ony)+.5d0)**2d0+
     +(abs(onz)+.5d0)**2d0)*l
      rr=0d0
      if(nx.ne.0) then
         rr=(abs(onx)-.5d0)**2d0
      endif
      if(ny.ne.0) then
         rr=rr+(abs(ony)-.5d0)**2d0
      endif
      if(nz.ne.0) then
         rr=rr+(abs(onz)-.5d0)**2d0
      endif
      rmin=sqrt(rr)*l
      return
      end
      logical function grid_incube(x,y,z,nx,ny,nz,l_rect)
      double precision x,y,z ! [in] point coors
      double precision l_rect ! [in] cube size
      integer nx,ny,nz       ! [in] cube definition
c
c     returns true if xyz is inside cube defined by nx,ny,nz
c
#include "cdft.fh"      
c     
      double precision onx,ony,onz

c
      onx=dble(nx)
      ony=dble(ny)
      onz=dble(nz)
      grid_incube=.false.
      if(x.le.(onx+.5d0)*l_rect.and.x.gt.(onx-.5d0)*l_rect) then
         if(y.le.(ony+.5d0)*l_rect.and.y.gt.(ony-.5d0)*l_rect) then
            if(z.le.(onz+.5d0)*l_rect.and.z.gt.(onz-.5d0)*l_rect) then
               grid_incube=.true.
            endif
         endif
      endif
      return
      end
      double precision function grid_lcube(rad)
      double precision rad
c
      double precision l_cube
      parameter(l_cube=4d0)
      integer n
c
      n=(rad+rad)/l_cube
      grid_lcube=dble(rad+rad)/dble(n)
      return
      end
