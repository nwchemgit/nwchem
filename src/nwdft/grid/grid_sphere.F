      logical function grid_insidesp(nx,ny,nz,l_rect,rad)
c
c$Id: grid_sphere.F,v 1.2 2000-02-23 01:37:57 edo Exp $
c
      implicit none
      integer nx,ny,nz
      double precision l_rect,rad
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
      integer function grid_ncubes(rad,ynxyz,fnxyz)
      implicit none
      integer nx,ny,nz,nxyz
      integer ninside
      double precision rad
      integer fnxyz(*)
      logical ynxyz
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
                     call n2nxyz(nxyz,fnxyz(ninside),onx,ony,onz)
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
      subroutine n2nxyz(nxyz,npoint,nx,ny,nz)
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
      integer function grid_wcube(x,y,z,fnxyz,rad,ncubes)
      integer ncubes
      double precision x,y,z, rad
      double precision fnxyz(*) ! remember to pass it as f(ntypes,1)
#include "cdft.fh"      
c     
      integer nx,ny,nz,nn
      integer lrec2
      integer nxyz2n,onx,ony,onz,onxyz
      nxyz2n(onx,ony,onz,onxyz)=(onx+onxyz)+
     +( (ony+onxyz)+(onz+onxyz)*(2*onxyz+1) )*(2*onxyz+1)
c
      nxyz=nint((rad*2d0-l_rect)/(l_rect*2d0))

c     get the cube
c
      lrec2=l_rec*2d0
      nx=nint(x/lrect2)
      ny=nint(y/lrect2)
      nz=nint(z/lrect2)
      nn=nxyz2n(nx,ny,nz,nxyz)
c
c     see which cube is associated by nn
c
      do i=1,ncubes
         if(fnxyz(i).eq.nn) then
            grid_wcube=i
            return
         endif
      enddo
      call errquit(' grid_wcube: failed',ncubes)
      return
      end
      
