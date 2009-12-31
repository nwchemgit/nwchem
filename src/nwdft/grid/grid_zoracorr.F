czora...
czora... Calculate the zora correction on the grid
czora...
      subroutine grid_zoracorr(nqpts,qxyz,qwght,natoms,g_dens,amat)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
c
      integer nqpts
      integer g_dens(2),igrid,natoms,npol
      double precision qxyz(*),qwght(*)![in]
      double precision amat(nqpts,ipol)
      double precision amat_coul(nqpts,ipol)
      double precision amat_nucl(nqpts)
      double precision tol
      double precision zoraCorr, totPot
      integer iatom
      double precision nucCharge, nucCoords(3)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision denomFac
      integer closegridpts(nqpts)
c
      double precision clight_au2
      clight_au2 = clight_au*clight_au
c
c     == preliminaries ==
      do igrid = 1,nqpts
        amat(igrid,1) = 0.d0            
        amat_coul(igrid,1) = 0.d0            
        if (ipol.gt.1) then
           amat(igrid,2) = 0.d0            
           amat_coul(igrid,2) = 0.d0            
        end if 
        amat_nucl(igrid) = 0.d0
        closegridpts(igrid) = 0
      end do
c
c     == calculate the hartree potential on a supplied list of points ==
      tol = 1d-8
      call potential_list(ao_bas_han, g_dens(1), nqpts,
     &     qxyz, amat_coul(1,1), tol)
    
      if (ipol.gt.1) then
        call potential_list(ao_bas_han, g_dens(2), nqpts,
     &     qxyz, amat_coul(1,2), tol)
      end if
c
c     == calculate the total nuclear potential on the grid ==
      call gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &    closegridpts,amat_nucl)
c
c     == assemble zora correction ==
      zoraCorr = 0.d0
      totPot = 0.d0
      do igrid = 1,nqpts
       if (ipol.gt.1) then
        totPot = - amat_coul(igrid,1) - amat_coul(igrid,2) 
     &         + amat_nucl(igrid)
       else
        totPot = - amat_coul(igrid,1) + amat_nucl(igrid)
       end if
c
c      == assemble the appropriate correction ==
       if (zora_calc_type.eq.0) then  ! pure kinetic test
          zoraCorr = 0.5d0
       else if (zora_calc_type.eq.1) then ! zora correction
             zoraCorr = totPot/(4.d0*clight_au2 - 2.d0*totPot)
       else if (zora_calc_type.eq.2) then ! zora energy scaling
             denomFac = (2.d0*clight_au2 - totPot)
             zoraCorr = clight_au2/denomFac/denomFac
       else
          zoraCorr = 0.d0
       end if
c
c      == multiply by the integration weight ==
       if (igrid.eq.closegridpts(igrid)) then
          zoraCorr = 0.d0
       else
          zoraCorr = zoraCorr*qwght(igrid)
       end if
c
c      == fill array ==
       amat(igrid,1) = zoraCorr
       if (ipol.gt.1) amat(igrid,2) = zoraCorr
c
      end do
c
      return
      end
c
czora...
czora...Calculate the nuclear potential on the grid
czora...
      subroutine gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &   closegridpts,amat_nucl)
c
      implicit none
c
#include "geom.fh"
#include "stdio.fh"
#include "zora.fh"
c
      integer iatom,natoms,geom,igrid,nqpts
      double precision nucCharge, nucCoords(3)
      double precision qxyz(3,nqpts),qwght(nqpts)
      character*16 tags(natoms)
      logical lSuccess
      double precision rx,ry,rz,dist
      double precision amat_nucl(nqpts)
      integer closegridpts(*)
c
c     write(luout,*) "zoracutoff:",zoracutoff
c
c     == get the total nuclear potential on the grid ==
c
c     == loop over the grid points ==
      do igrid = 1,nqpts
c
      amat_nucl(igrid) = 0.d0
c     == loop over the atoms ==
      do iatom = 1,natoms
c
c     == get an atom (needs error handling) ==
      lSuccess = geom_cent_get(geom, iatom, tags, nucCoords, 
     &          nucCharge)
c
c     == distance from the grid points to the atom centers == 
      rx = nucCoords(1) - qxyz(1,igrid)
      ry = nucCoords(2) - qxyz(2,igrid)
      rz = nucCoords(3) - qxyz(3,igrid)
      dist = dsqrt(rx*rx + ry*ry + rz*rz)
c
c     == check cutoff ==
      if (dist.gt.zoracutoff) then
         amat_nucl(igrid) = amat_nucl(igrid) - nucCharge/dist
      else
         closegridpts(igrid) = igrid
      end if
c
      end do
      end do
c
      return
      end 
