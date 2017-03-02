c
      subroutine calc_zora_facs(ao_bas_han,geom,
     &                          ipol,g_dens,natoms,npts,qxyz,qwght,
     &                          use_modelpotential,gexpo,gcoef,
     &                          ofinite,    ! in: =.true. if finite Nuclear Model (Gaussian)
     &                          zetanuc_arr,! in: zetanuc(i) i=1,natoms
     &                          Knucl,      ! in: = .true. if including ONLY nuclear part in K ZORA
     &                          fac1_arr,
     &                          fac2_arr)
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
      double precision qwght(npts)
      double precision pot(npts)
      double precision qxyz(3,npts)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts),fac2_arr(npts)
      logical use_modelpotential
      double precision gexpo(natoms,50)
      double precision gcoef(natoms,50)
      logical ofinite,Knucl
      double precision zetanuc_arr(natoms)
c
c     == preliminaries ==
      clight_au2 = clight_au*clight_au
      do ipt = 1,npts
        amat_coul(ipt,1) = 0.d0            
        if (ipol.gt.1) amat_coul(ipt,2) = 0.d0            
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do
c
c     == calculate the total hartree potential on the grid ==
      call gridHartreePotential(use_modelpotential,
     &    ao_bas_han, geom, natoms, ipol, g_dens, npts, qxyz, qwght,
     &    closegridpts, gexpo, gcoef, amat_coul)
c
c     == calculate the total nuclear potential on the grid ==
      if (ofinite) then
c         == Choosing Nuclear Model: erf(zetanuc^0.5 r_L) ==
          call gridNuclearPotentialFinite(geom,natoms,npts,qxyz,qwght,
     &                             zetanuc_arr,
     &                             closegridpts,amat_nucl)
      else ! default : point charge model for nuclei
          call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                             closegridpts,amat_nucl)
      end if !ofinite
c
c     == zero out weights of close points ==
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c
c     == assemble zora factors ==
      if (Knucl) then ! V=Vnucl
       do k = 1,npts
         amat_tot = amat_nucl(k) ! pure nuclear
         fac1_arr(k)=amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
     &              *qwght(k)
         denom = (2.d0*clight_au2 - amat_tot)
         fac2_arr(k)=(clight_au2/denom/denom)
     &              *qwght(k)
       enddo ! end loop k grid points
      else            ! V=Vnucl+Vee (default)
       do k = 1,npts
c       == assemble hartree and nuclear contributions ==
        amat_tot = amat_coul(k,1)+amat_nucl(k)
        fac1_arr(k)=amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
     &              *qwght(k)
        denom = (2.d0*clight_au2 - amat_tot)
        fac2_arr(k)=(clight_au2/denom/denom)
     &              *qwght(k)
       enddo ! end loop k grid points
      endif ! end-if-Knucl
c
      return
      end
