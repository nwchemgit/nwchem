      logical function grid_insidesp(fnxyz,l_rect,rad)
c
c$Id: grid_sphere.F,v 1.11 2000-04-06 01:25:44 edo Exp $
c
      implicit none
      integer fnxyz ! [in]
      double precision l_rect,rad ! [in]
c
      integer nx,ny,nz            
      double precision ndx,ndy,ndz,dista
      integer nxyz,nxyz2
      integer grid_nxyz
      double precision ll
      external grid_nxyz
c
#include "grid_cube.fh"      
c
      grid_insidesp=.true.
      nxyz=grid_nxyz(rad,l_rect)
      call grid_n2nxyz(nxyz,fnxyz,nx,ny,nz)
      if(fnxyz.lt.shft.and.sub4) then
         if(nx.eq.2) then
            ndx=0d0
         else
            ndx=0.5d0
         endif
         if(ny.eq.2) then
            ndy=0d0
         else
            ndy=0.5d0
         endif
         if(nz.eq.2) then
            ndz=0d0
         else
            ndz=0.5d0
         endif
         ll=l_rect*2d0/3d0
      else
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
         ll=l_rect
      endif
      dista=ll*sqrt(ndx*ndx+ndy*ndy+ndz*ndz)
      if(dista.gt.rad) grid_insidesp=.false.
      return
      end
      integer function grid_ncubes(rad,ynxyz,fnxyz,l_rect)
      implicit none
#include "cgridfile.fh"
#include "cdft.fh"
      integer nx,ny,nz,nxyz
      integer ninside
      double precision rad ! [in]
      integer fnxyz(*)     ! [out]
      logical ynxyz        ! [in] get points
      logical grid_insidesp
      double precision l_rect ! [in]
      external grid_insidesp
      integer grid_nxyz,ffnxyz,ncount
      external grid_nxyz
      integer nxyz2n,onx,ony,onz,onxyz,nxyz2,ii,nchk
c
#include "grid_cube.fh"      
      nxyz2n(onx,ony,onz,onxyz,shft)=onx-1+
     +(ony-1+(onz-1)*onxyz)*onxyz+shft
c
      nxyz=grid_nxyz(rad,l_rect)
      nxyz2=nxyz/2
      ninside=0
      if(sub4) then
      ncount=-1
c     inside part of side L/2
      do nz=1,4
         do ny=1,4
            do nx=1,4
               ncount=ncount+1
               ffnxyz=ncount
               if(grid_insidesp(ffnxyz,l_rect,rad)) then
                  ninside=ninside+1
                  if(ynxyz) then
                     fnxyz(ninside)=ncount
!               write(0,'("A0",2I4)') ninside,fnxyz(ninside)
                     call grid_n2nxyz(nxyz,fnxyz(ninside),onx,ony,onz)
                     if(onx.ne.nx.or.ony.ne.ny.or.onz.ne.nz) then
                        write(0,'(A,4I4)') ' IN  nx ny nz ',nx,ny,nz,
     ,                       fnxyz(ninside)
                        write(0,'(A,4I4)') ' OUT nx ny nz ',onx,ony,onz
                        call errquit('grid_ncube: mapping wrong',9)
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
      endif
      do nx=1,nxyz
         do ny=1,nxyz
            do nz=1,nxyz
               if(sub4) then
                  if((nx.eq.nxyz2.or.nx.eq.nxyz2+1).and.
     .                 (ny.eq.nxyz2.or.ny.eq.nxyz2+1).and.
     .                 (nz.eq.nxyz2.or.nz.eq.nxyz2+1)) goto 413
c     do nothing --- see above
               endif
               ffnxyz=nxyz2n(nx,ny,nz,nxyz,shft)
               if(grid_insidesp(ffnxyz,l_rect,rad)) then
                  ninside=ninside+1
                  if(ynxyz) then
                     fnxyz(ninside)=ffnxyz
c     write(0,'("A1",2I4)') ninside,fnxyz(ninside)
                     nchk=fnxyz(ninside)
                     call grid_n2nxyz(nxyz,nchk,onx,ony,onz)
                     if(onx.ne.nx.or.ony.ne.ny.or.onz.ne.nz) then
                        write(0,'(A,4I4)') ' IN  nx ny nz ',nx,ny,nz,
     ,                       fnxyz(ninside)
                        write(0,'(A,4I4)') ' OUT nx ny nz ',onx,ony,onz
                        call errquit('grid_ncube: mapping wrong',9)
                     endif
                  endif
               endif
 413           continue
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
      integer nresid,factxyz,npu
#include "grid_cube.fh"      
c      
      if(sub4.and.npoint.lt.shft) then
         npu=npoint
         factxyz=4
      else
         factxyz=nxyz
!         npu=npoint
         npu=npoint-shft
      endif
      nz=npu/(factxyz*factxyz)+1
      nresid=npu-(nz-1)*factxyz*factxyz
      ny=nresid/factxyz+1
      nresid=nresid-(ny-1)*factxyz
      nx=nresid+1
      return
      end
      subroutine grid_rmaxmin(fnxyz,rmax,rmin,l,nxyz)
      implicit none
      integer fnxyz      ! [in]
      integer nxyz ! [in]
      double precision l    ! [in]
      double precision rmax ! [out]
      double precision rmin ! [out]
c
      double precision onx,ony,onz,nxyz2
      double precision onxyz2
      integer nx,ny,nz      
c
#include "grid_cube.fh"      
c
      rmax=0d0
      rmin=0d0
      call grid_n2nxyz(nxyz,fnxyz,nx,ny,nz)
      onx=nx
      ony=ny
      onz=nz
      if(sub4.and.fnxyz.lt.shft) then
         if(nx.le.2) then
            rmax=rmax+(3d0-onx)**2
            rmin=rmin+(2d0-onx)**2
         else
            rmax=rmax+(onx-2d0)**2
            rmin=rmin+(onx-3d0)**2
         endif
         if(ny.le.2) then
            rmax=rmax+(3d0-ony)**2
            rmin=rmin+(2d0-ony)**2
         else
            rmax=rmax+(ony-2d0)**2
            rmin=rmin+(ony-3d0)**2
         endif
         if(nz.le.2) then
            rmax=rmax+(3d0-onz)**2
            rmin=rmin+(2d0-onz)**2
         else
            rmax=rmax+(onz-2d0)**2
            rmin=rmin+(onz-3d0)**2
         endif
         rmax=sqrt(rmax)*l*0.5d0
         rmin=sqrt(rmin)*l*0.5d0
      else
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
      endif
      return
      end
      logical function grid_incube(x,y,z,nx,ny,nz,fnxyz,nxyz2,ll)
      double precision x,y,z ! [in] point coors
      integer nx,ny,nz       ! [in] cube definition
      double precision ll ! [in]
      integer fnxyz ! [in]
      integer nxyz2 ! [in]
c
c     returns true if xyz is inside cube defined by nx,ny,nz
c
#include "cdft.fh"      
#include "grid_cube.fh"      
c     
      double precision l23
c
      grid_incube=.false.
      if(sub4.and.fnxyz.lt.shft) then
         l23=ll*0.5d0
         if(x.gt.(nx-3d0)*l23.and.x.le.(nx-2d0)*l23) then
            if(y.gt.(ny-3d0)*l23.and.y.le.(ny-2d0)*l23) then
               if(z.gt.(nz-3d0)*l23.and.z.le.(nz-2d0)*l23) then
                  grid_incube=.true.
               endif
            endif
         endif
      else
         if(x.le.(nx-nxyz2)*ll.and.x.gt.(nx-nxyz2-1)*ll) then
            if(y.le.(ny-nxyz2)*ll.and.y.gt.(ny-nxyz2-1)*ll) then
               if(z.le.(nz-nxyz2)*ll.and.z.gt.(nz-nxyz2-1)*ll) then
                  grid_incube=.true.
               endif
            endif
         endif
      endif
!      write(0,'(A,5I5,3F12.6)') ' fnxyz nx ny nz ',
!     ,     fnxyz,nx,ny,nz,grid_incube,
!     ,     x,y,z
c      write(0,1) x,y,z,nx,ny,nz,grid_incube
c 1    format(3f20.12,3I3,I5)
      return
      end
      double precision function grid_lcube(rad)
      implicit none
#include "grid_cube.fh"      
      double precision rad
c
c
      integer n
c
      if(l_cube.lt.0d0) then
         n=-int(l_cube)
         grid_lcube=-dble(rad+rad)/l_cube
      else
         n=max(nint((rad+rad)/l_cube),2)
         grid_lcube=dble(rad+rad)/dble(n)
      endif
      return
      end
      integer function grid_nxyz(rad,l_rect)
      implicit none
      double precision rad,l_rect
      grid_nxyz=nint(rad/l_rect)*2d0
      return
      end
