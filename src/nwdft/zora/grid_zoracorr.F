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
#include "global.fh"  ! FA
c
      external gridQpqPotential
      integer nqpts
      integer g_dens(2),igrid,natoms,npol
      double precision qxyz(3,nqpts),qwght(nqpts) ![in]
      double precision amat(nqpts,ipol)
      double precision amat_coul(nqpts,ipol)
      double precision amat_nucl(nqpts)
      double precision amat_Qnucl(nqpts) ! Added by FA
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
        amat_Qnucl(igrid) = 0.d0 ! Added by FA
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
     &                          closegridpts,amat_nucl)
      if ((zora_calc_type.eq.3).or.(zora_calc_type.eq.4)) then
c     == calculate Quadrupole potential on the grid ==
       call gridQpqPotential(nqpts,qxyz,amat_Qnucl,
     &                       closegridpts)   
      end if
c
c     == assemble zora correction ==
      zoraCorr = 0.d0
      totPot = 0.d0
      do igrid = 1,nqpts
       if (ipol.gt.1) then
        totPot = -amat_coul(igrid,1)-amat_coul(igrid,2) 
     &           + amat_nucl(igrid)
       else
        totPot = -amat_coul(igrid,1)+amat_nucl(igrid)
       end if
c
c      == assemble the appropriate correction ==
       if      (zora_calc_type.eq.0) then  ! pure kinetic test
             zoraCorr = 0.5d0
       else if (zora_calc_type.eq.1) then ! zora correction
             zoraCorr = totPot/(4.d0*clight_au2 - 2.d0*totPot)
       else if (zora_calc_type.eq.2) then ! zora energy scaling
             denomFac = (2.d0*clight_au2 - totPot)
             zoraCorr = clight_au2/denomFac/denomFac
       else if (zora_calc_type.eq.3) then ! zora EFG      - Added by FA
             denomFac = (2.d0*clight_au2 - totPot)
             zoraCorr = clight_au2/denomFac/denomFac
             zoraCorr = zoraCorr*amat_Qnucl(igrid)
       else if (zora_calc_type.eq.4) then ! Numerical EFG - Added by FA
             zoraCorr = amat_Qnucl(igrid)
c            zoraCorr = 1.0   ! For Atomic Nr FA-04-05-10
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
      end do
c
      return
      end
c $Id$
