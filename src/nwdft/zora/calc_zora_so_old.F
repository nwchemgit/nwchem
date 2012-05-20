c
c     == calculate spin-orbit zora contribution ==
      subroutine calc_zora_so(ao_bas_han, geom, ipol, g_dens, 
     &      chi_ao, delchi_ao, qxyz, qwght, nbf, npts, natoms, 
     &      ovlp, ke, 
     &      zora0, scal0, 
     &      zorax, scalx, 
     &      zoray, scaly, 
     &      zoraz, scalz)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"
c
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision chi_ao(npts,nbf)
      double precision qwght(npts)
      double precision pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
c
      double precision ovlp(nbf,nbf)
      double precision ke(nbf,nbf)
      double precision zora0(nbf,nbf)      ! scalar component
      double precision scal0(nbf,nbf)  ! scalar component
      double precision zorax(nbf,nbf)
      double precision scalx(nbf,nbf)
      double precision zoray(nbf,nbf)
      double precision scaly(nbf,nbf)
      double precision zoraz(nbf,nbf)
      double precision scalz(nbf,nbf)
c
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot,fac1,fac2
      double precision denom
c
      clight_au2 = clight_au*clight_au
c
c     == preliminaries ==
      do ipt = 1,npts
        amat_coul(ipt,1) = 0.d0            
        if (ipol.gt.1) amat_coul(ipt,2) = 0.d0            
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do
c
c     == calculate the hartree potential on a supplied list of points ==
      tol = 1d-8
      call potential_list(ao_bas_han, g_dens(1), npts, qxyz, 
     &  amat_coul(1,1), tol)
      if (ipol.gt.1) then
        call potential_list(ao_bas_han, g_dens(2), npts, qxyz, 
     &  amat_coul(1,2), tol)
      end if
c
c     == calculate the total nuclear potential on the grid ==
      call gridNuclearPotential(geom,natoms,npts,qxyz,qwght,
     &                          closegridpts,amat_nucl)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c
c     == assemble zora correction ==
      amat_tot = 0.d0
      do i = 1, nbf
        do j = 1, nbf
          do k = 1, npts
c
c          == assemble hartree and nuclear contributions ==
           if (ipol.gt.1) then
             amat_tot = -amat_coul(k,1)-amat_coul(k,2)+amat_nucl(k)
           else
             amat_tot = -amat_coul(k,1)+amat_nucl(k)
           end if
c
c          == calculate the overlap (for testing) ==
c          ovlp(i,j) = ovlp(i,j)+chi_ao(k,i)*qwght(k)*chi_ao(k,j)
c
c          == calculate the pure kinetic contribution (for testing) ==
c          ke(i,j) = ke(i,j)
c     &         +qwght(k)*(delchi_ao(k,1,i)*delchi_ao(k,1,j)
c     &                   +delchi_ao(k,2,i)*delchi_ao(k,2,j)
c     &                   +delchi_ao(k,3,i)*delchi_ao(k,3,j))
c
c          == calculate scalar zora and spin-orbit corrections ==
           fac1 = amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
           zora0(i,j) = zora0(i,j)
     &         +qwght(k)*(delchi_ao(k,1,i)*fac1*delchi_ao(k,1,j)
     &                   +delchi_ao(k,2,i)*fac1*delchi_ao(k,2,j) 
     &                   +delchi_ao(k,3,i)*fac1*delchi_ao(k,3,j))   
           zorax(i,j) = zorax(i,j)
     &         +qwght(k)*(delchi_ao(k,2,i)*fac1*delchi_ao(k,3,j)
     &                   -delchi_ao(k,3,i)*fac1*delchi_ao(k,2,j)) 
           zoray(i,j) = zoray(i,j)
     &         +qwght(k)*(delchi_ao(k,3,i)*fac1*delchi_ao(k,1,j)
     &                   -delchi_ao(k,1,i)*fac1*delchi_ao(k,3,j)) 
           zoraz(i,j) = zoraz(i,j)
     &         +qwght(k)*(delchi_ao(k,1,i)*fac1*delchi_ao(k,2,j)
     &                   -delchi_ao(k,2,i)*fac1*delchi_ao(k,1,j)) 
c
c          == calculate scalar zora and spin-orbit scaling ==
           denom = (2.d0*clight_au2 - amat_tot)
           fac2 =  clight_au2/denom/denom
           scal0(i,j) = scal0(i,j)
     &         +qwght(k)*(delchi_ao(k,1,i)*fac2*delchi_ao(k,1,j)
     &                   +delchi_ao(k,2,i)*fac2*delchi_ao(k,2,j) 
     &                   +delchi_ao(k,3,i)*fac2*delchi_ao(k,3,j))   
           scalx(i,j) = scalx(i,j)
     &         +qwght(k)*(delchi_ao(k,2,i)*fac2*delchi_ao(k,3,j)
     &                   -delchi_ao(k,3,i)*fac2*delchi_ao(k,2,j)) 
           scaly(i,j) = scaly(i,j)
     &         +qwght(k)*(delchi_ao(k,3,i)*fac2*delchi_ao(k,1,j)
     &                   -delchi_ao(k,1,i)*fac2*delchi_ao(k,3,j)) 
           scalz(i,j) = scalz(i,j)
     &         +qwght(k)*(delchi_ao(k,1,i)*fac2*delchi_ao(k,2,j)
     &                   -delchi_ao(k,2,i)*fac2*delchi_ao(k,1,j)) 
c
          enddo
        enddo
      enddo
c
      return
      end
c $Id$
