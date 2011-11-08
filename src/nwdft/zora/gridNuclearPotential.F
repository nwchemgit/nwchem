c
c     == Calculate the point nuclear potential on the grid ==
      subroutine gridNuclearPotentialPoint(geom,
     &                                natoms,      ! IN : number of atoms
     &                                nqpts,       ! IN : number of grid points
     &                                qxyz,        ! IN : grid points
     &                                qwght,       ! IN : quadrature weightings
     &                                closegridpts,
     &                                amat_nucl)   ! OUT : nuclear potential
      implicit none
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh" ! Added for using ga_nodeid()

      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts)
      integer closegridpts(*)
      logical ofinite
c     == get the total nuclear potential on the grid ==
      do igrid = 1,nqpts   ! == loop over the grid points ==
       amat_nucl(igrid) = 0.d0  
       do iatom = 1,natoms ! == loop over the atoms ==
c
c     == get an atom (needs error handling) ==
        lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &             nucCharge)
c
c     == distance from the grid points to the atom centers == 
        rx = nucCoords(1) - qxyz(1,igrid)
        ry = nucCoords(2) - qxyz(2,igrid)
        rz = nucCoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff) then ! == check cutoff ==
           amat_nucl(igrid) = amat_nucl(igrid) - 
     &                        nucCharge/dist
        else
          closegridpts(igrid) = igrid
        end if
       end do
      end do
      return
      end 
c
c     == Finite nucleus 1 ==
      subroutine gridNuclearPotentialFinite(
     &               geom,
     &               natoms,      ! IN : number of atoms
     &               nqpts,       ! IN : number of grid points
     &               qxyz,        ! IN : grid points
     &               qwght,       ! IN : quadrature weightings
     &               zetanuc_arr, ! IN : sqrt(zetanuc) for Gaussian Nuclear Model
     &               closegridpts,
     &               amat_nucl)   ! OUT : nuclear potential
c
c Evaluates: V_L^{Gauss}(r) = \sum_L -Z_L/r_L erf(\tilde{r}_L)
c            \tilde{r}_L = srqt(zetanuc) r_L
c            r_L = |r_grid-r_N|
c            r_N   , nuclear position
c            r_grid, grid point
c            according to: Autschbach, J. ChemPhysChem 2009,V10,P 22774-2283
c            page 2276
c Author: Fredy W. Aquino 
c Date  : 04-20-11
c
      implicit none
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh" ! Added for using ga_nodeid()

      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts),
     &                 rtemp
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts),
     &                zetanuc_arr(natoms),erf
      integer closegridpts(*)
      logical ofinite
      external merf,get_zetanuc_arr
      do igrid = 1,nqpts   ! == loop over the grid points ==
       amat_nucl(igrid) = 0.d0  
       do iatom = 1,natoms ! == loop over the atoms ==
c     == get an atom (needs error handling) ==
        lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &             nucCharge)
c     == distance from the grid points to the atom centers == 
        rx = nucCoords(1) - qxyz(1,igrid)
        ry = nucCoords(2) - qxyz(2,igrid)
        rz = nucCoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        rtemp = zetanuc_arr(iatom)*dist
        call merf(rtemp,erf) ! Evaluate erf(sqrt(zetanuc)r_L)
        if (dist.gt.zoracutoff) then ! == check cutoff ==          
           amat_nucl(igrid) = amat_nucl(igrid) - 
     &                        nucCharge/dist*erf
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-loop-iatom
      end do ! end-loop-igrid
      return
      end
c
c     == Finite nucleus 2 ==
      subroutine gridNuclearPotentialFinite2(
     &               geom,
     &               natoms,      ! IN : number of atoms
     &               nqpts,       ! IN : number of grid points
     &               qxyz,        ! IN : grid points
     &               qwght,       ! IN : quadrature weightings
     &               closegridpts,
     &               amat_nucl)   ! OUT : nuclear potential
c Evaluates: V_L^{Gauss}(r) = - Z_L/r_L P(1/2,\tilde{r}_L^2)
c            \tilde{r}_L = srqt(zetanuc) r_L
c            r_L = |r_grid-r_N|
c            r_N   , nuclear position
c            r_grid, grid point
c            according to: Autschbach, J. ChemPhysChem 2009,V10,P 22774-2283
c            Eq. (5). 
c Author: Fredy W. Aquino 
c Date  : 04-20-11
      implicit none
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts),
     &                zetanuc_arr(natoms)
      real a_coeff,rtemp,gammap,gammaq
      integer closegridpts(*)
      external gratio,  ! Evaluates Incomplete Gamma Function P(a,x)
     &         get_zetanuc_arr

      call get_zetanuc_arr(geom,natoms,zetanuc_arr) ! zetanuc_arr(i) i=1,natoms

      a_coeff=0.5
      do igrid = 1,nqpts   ! == loop over the grid points ==
       amat_nucl(igrid) = 0.d0  
       do iatom = 1,natoms ! == loop over the atoms ==
c     == get an atom (needs error handling) ==
        lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &             nucCharge)
c     == distance from the grid points to the atom centers == 
        rx = nucCoords(1) - qxyz(1,igrid)
        ry = nucCoords(2) - qxyz(2,igrid)
        rz = nucCoords(3) - qxyz(3,igrid)
        dist = dsqrt(rx*rx + ry*ry + rz*rz)
        rtemp = zetanuc_arr(iatom)*dist*dist
        call gratio(a_coeff,rtemp,gammap,gammaq,0)! Evaluate P(1/2,zetanuc r_L^2)
        if (dist.gt.zoracutoff) then ! == check cutoff ==          
           amat_nucl(igrid) = amat_nucl(igrid) - 
     &                        nucCharge/dist*gammap
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-loop-iatom
      end do ! end-loop-igrid
      return
      end
c $Id$
