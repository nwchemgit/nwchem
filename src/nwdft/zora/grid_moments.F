 
      subroutine grid_moments(
     &   origin,                ! in : origin
     &   type,                  ! in: integral type
     &   ncomp,                 ! in: no of components
     &   nqpts,                 ! in : nr. grid points 
     &   qxyz,                  ! in : grid points             
     &   amat,                  ! out: moments on grid
     &   closegridpts)

c           nqpts    , number of grid points
c           qxyz     , grid points

      implicit none
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"
#include "errquit.fh"
      integer igrid,nqpts,ncomp,dbg
      integer closegridpts(*)
      double precision origin(3)
      double precision qxyz(3,nqpts)
      character*(*) type
      double precision rx,ry,rz,dist
      double precision amat(nqpts,ncomp) ! output
      logical master

c     =================================================================

      dbg = 0
      master = (ga_nodeid().eq.0)
      if (dbg>0 .and. master) write (luout,*) 'grid_moments: cutoff =',
     &   zoracutoff_EFG


c     == loop over the grid points ==
      
      do igrid = 1,nqpts
c       == distance from the grid points to given origin == 
        rx = qxyz(1,igrid)-origin(1)
        ry = qxyz(2,igrid)-origin(2)
        rz = qxyz(3,igrid)-origin(3)
c        dist = dsqrt(rx*rx + ry*ry + rz*rz)

        if (ncomp.eq.3) then
c         dipole
          amat(igrid,1) = rx
          amat(igrid,2) = ry
          amat(igrid,3) = rz
        else
          call errquit("grid_moments: ncomp unsupported", 666,
     &       INPUT_ERR)
        end if

      end do                    ! end-grid
      
      return
      end 
