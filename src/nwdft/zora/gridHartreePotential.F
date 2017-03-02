c
c     Calculate the Hartree potential on the grid
      subroutine gridHartreePotential(use_modelpotential,
     &       ao_bas_han, geom, natoms, ipol, g_dens, npts, qxyz, qwght,
     &       closegridpts, gexpo, gcoef, amat_coul)
c
      implicit none
c
#include "stdio.fh"
c
      integer npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision qxyz(3,npts)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      integer closegridpts(npts)
      double precision tol
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
c
c     == use the model potential approach ==
      if (use_modelpotential) then
c
c      == calculate the model density to construct the Coulomb contribution ==
       call calc_modelpotential(geom,natoms,npts,qxyz,qwght,
     &       closegridpts,gexpo,gcoef,amat_coul)  !amat_coul = v(modelpot)
      else
c
c      == calculate the hartree potential on a supplied list of points ==
       tol = 1d-8
       do i=1,ipol
          call potential_list(ao_bas_han, g_dens(i), npts, qxyz,
     &                     amat_coul(1,i), tol)
       enddo
       do k = 1,npts
          if (ipol.gt.1) amat_coul(k,1)=amat_coul(k,1)+amat_coul(k,2)
          amat_coul(k,1) = -amat_coul(k,1)
       enddo
      end if ! use_modelpotential

      return
      end
