      subroutine gridQpqPotential(zora_Qpq,      ! in : type of EFG potential
     &                            xyz_EFGcoords, ! in : EFG-nuclear coord.
     &                            nqpts,         ! in : nr. grid points 
     &                            qxyz,          ! in : grid points             
     &                            amat_Qnucl,    ! out: EFG potential
     &                            closegridpts)
c    Purpose: Evaluates Q_pq(\vec{r},\vec{r}')
c    Input: 
c           zora_Qpq = 1  -> Evaluates Qxx
c                    = 2  -> Evaluates Qyy
c                    = 3  -> Evaluates Qzz
c                    = 4  -> Evaluates Qxy
c                    = 5  -> Evaluates Qxz
c                    = 6  -> Evaluates Qyz
c           nqpts    , number of grid points
c           qxyz     , grid points
c           xyz_EFGcoords, Quadrupole potential is evaluated
c                          in this point (\vec{r})
c    Output:
c           amat_Qnucl, Quadrupole potential ev. in the grid
c                       for integration purpose (\vec{r}')
c    Author: Fredy Aquino
c    Date  : 10-29-10

      implicit none
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"
      integer igrid,nqpts
      integer closegridpts(*)
      integer zora_Qpq
      double precision xyz_EFGcoords(3)
      double precision qxyz(3,nqpts)
      double precision rx,ry,rz,dist,dist5
      double precision amat_Qnucl(nqpts)

c     == loop over the grid points ==
      if (zora_Qpq.eq.1) then        ! Qxx
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
          amat_Qnucl(igrid) = (2.d0*rx*rx-(ry*ry+rz*rz))/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      else if (zora_Qpq.eq.2) then    ! Qyy
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
         amat_Qnucl(igrid) = (2.d0*ry*ry-(rx*rx+rz*rz))/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      else if (zora_Qpq.eq.3) then    ! Qzz
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
         amat_Qnucl(igrid) = (2.d0*rz*rz-(rx*rx+ry*ry))/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      else if (zora_Qpq.eq.4) then    ! Qxy
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
          amat_Qnucl(igrid) = (3.d0*rx*ry)/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      else if (zora_Qpq.eq.5) then    ! Qxz
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
          amat_Qnucl(igrid) = (3.d0*rx*rz)/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      else if (zora_Qpq.eq.6) then    ! Qyz
       do igrid = 1,nqpts
c     == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then ! check-cutoff
          dist5=dist*dist*dist*dist*dist
          amat_Qnucl(igrid) = (3.d0*ry*rz)/dist5
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      endif
      return
      end 


      subroutine gridQpqPotential_all(
     &   xyz_EFGcoords,         ! in : EFG-nuclear coord.
     &   nqpts,                 ! in : nr. grid points 
     &   qxyz,                  ! in : grid points             
     &   amat_Qnucl,            ! out: EFG potential
     &   closegridpts)
c    Purpose: Evaluates Q_pq(\vec{r},\vec{r}')
c     for all 6 unique pq simultaneously
c             index  = 1  -> Evaluates Qxx
c                    = 2  -> Evaluates Qyy
c                    = 3  -> Evaluates Qzz
c                    = 4  -> Evaluates Qxy
c                    = 5  -> Evaluates Qxz
c                    = 6  -> Evaluates Qyz
c           nqpts    , number of grid points
c           qxyz     , grid points
c           xyz_EFGcoords, Quadrupole potential is evaluated
c                          in this point (\vec{r})
c    Output:
c           amat_Qnucl, Quadrupole potential ev. in the grid
c                       for integration purpose (\vec{r}')
c
c     Note: these are the traceless versions where xx is actually
c     3xx - rr, etc


      implicit none
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"
      integer igrid,nqpts
      integer closegridpts(*)
      double precision xyz_EFGcoords(3)
      double precision qxyz(3,nqpts)
      double precision rx,ry,rz,dist,dist5
      double precision amat_Qnucl(nqpts,6) ! output
      integer dbg
      logical me

c     =================================================================

      dbg = 0
      me = (ga_nodeid().eq.0)
      if (dbg>0 .and. me) write (luout,*) 'gridQpq: cutoff =',
     &   zoracutoff_EFG


c     == loop over the grid points ==
      
      do igrid = 1,nqpts
c       == distance from the grid points to given xyz_EFGcoords() == 
        rx = xyz_EFGcoords(1) - qxyz(1,igrid)
        ry = xyz_EFGcoords(2) - qxyz(2,igrid)
        rz = xyz_EFGcoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)

        if (dist.gt.zoracutoff_EFG) then ! check-cutoff

          dist5=dist*dist*dist*dist*dist
          
          amat_Qnucl(igrid,1) = (2.d0*rx*rx-(ry*ry+rz*rz))/dist5
          amat_Qnucl(igrid,2) = (2.d0*ry*ry-(rx*rx+rz*rz))/dist5
          amat_Qnucl(igrid,3) = (2.d0*rz*rz-(rx*rx+ry*ry))/dist5
          amat_Qnucl(igrid,4) = (3.d0*rx*ry)/dist5
          amat_Qnucl(igrid,5) = (3.d0*rx*rz)/dist5
          amat_Qnucl(igrid,6) = (3.d0*ry*rz)/dist5

        else
          closegridpts(igrid) = igrid
        end if
      end do                    ! end-grid
      
      return
      end 

c $Id: gridQpqPotential.F 22517 2012-05-24 18:37:27Z jochen $
