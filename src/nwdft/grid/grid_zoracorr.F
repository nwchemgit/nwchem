czora...
czora... Calculate the zora correction on the grid
czora...
      subroutine grid_zoracorr(nqpts,qxyz,qwght,natoms,g_dens,amat)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"

      integer nqpts
      integer g_dens(2),igrid,natoms
      double precision qxyz(*),qwght(*)![in]
      double precision amat(nqpts,ipol)
      double precision amat_coul(nqpts)
      double precision amat_nucl(nqpts)
      double precision tol
      double precision zoraCorr, totPot
      integer iatom
      double precision nucCharge, nucCoords(3)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist

      double precision clight_au2
      clight_au2 = clight_au*clight_au

c     initialisation
      call dfill(nqpts, 0.d0, amat_coul, 1)
      call dfill(nqpts, 0.d0, amat_nucl, 1)

c     calculate the hartree potential on a supplied list of points
      tol = 1d-8
      call potential_list(ao_bas_han, g_dens(1), nqpts,
     &     qxyz, amat_coul, tol)

c     calculate the total nuclear potential on the grid
      call gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &    amat_nucl)

c     assemble zora correction
      zoraCorr = 0.d0
      totPot = 0.d0
      do igrid = 1,nqpts
        totPot = -amat_coul(igrid) + amat_nucl(igrid)
        zoraCorr = totPot/(4.d0*clight_au2 - 2.d0*totPot)
        zoraCorr = zoraCorr*qwght(igrid)
c                                           zoraCorr = 0.5d0*qwght(igrid)  ! pure kinetic test
        amat(igrid,1) = zoraCorr
        if (ipol .gt. 1) amat(igrid,2) = zoraCorr
      end do

      return
      end

czora...
czora...Calculate the nuclear potential on the grid
czora...
      subroutine gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &    amat_nucl)

      implicit none

#include "geom.fh"

      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts)

c     get the total nuclear potential on the grid

c     loop over the atoms
      do iatom = 1,natoms

c       get an atom (needs error handling)
        lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &          nucCharge)

c       loop over the grid points
        do igrid = 1,nqpts

c           distance from the grid points to the atom centers       
            rx = nucCoords(1) - qxyz(1,igrid)
            ry = nucCoords(2) - qxyz(2,igrid)
            rz = nucCoords(3) - qxyz(3,igrid)
            dist = dsqrt(rx*rx + ry*ry + rz*rz)

c           needs a zero or small distance check here
            amat_nucl(igrid) = amat_nucl(igrid) 
     &                  - nucCharge/dist

        end do
      end do

      return
      end 
