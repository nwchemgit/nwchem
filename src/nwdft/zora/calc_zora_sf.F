c
c     == calculate spin-free zora contribution ==
      subroutine calc_zora_sf(ao_bas_han, ! in: AO basis handle 
     &                           geom,       ! in: geometry handle 
     &                           ipol,       ! in: nr. of polarizations 
     &                           g_dens,     ! in: superposit. atomic densities 
     &                           chi_ao,     ! in: basis functions
     &                           delchi_ao,  ! in: deriv. of basis functions
     &                           qxyz,       ! in: grid points 
     &                           qwght,      ! in: weighting coeffs.
     &                           nbf,        ! in: nr. basis functions
     &                           npts,       ! in: nr. grid points
     &                           natoms,     ! in: nr. atoms 
     &                           zorasf,     ! out
     &                           zorascal,   ! out
     &                           ofinite,    ! in: =.true. if finite Nuclear Model (Gaussian)
     &                           zetanuc_arr,! in: zetanuc(i) i=1,natoms
     &                           Knucl)      ! in: = .true. if including ONLY nuclear part in K ZORA  
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"

      integer nbf, npts, ao_bas_han, natoms, geom
      integer g_dens(2),ipol
      double precision chi_ao(npts,nbf),qwght(npts),pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      double precision zorasf(nbf,nbf)
      double precision zorascal(nbf,nbf)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl,prod
      logical ofinite,Knucl
c
c ------- for Gaussian Nuclear Model --- START
      double precision zetanuc_arr(natoms)
c ------- for Gaussian Nuclear Model --- START
      external get_ints_zora_sf,gridNuclearPotentialPoint,
     &         gridNuclearPotentialFinite,
     &         gridNuclearPotentialFinite2
c
c     == preliminaries ==
      clight_au2 = clight_au*clight_au
      do ipt = 1,npts
        do i=1,ipol
         amat_coul(ipt,i) = 0.d0            
        end do
        amat_nucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do
c
c     == calculate the hartree potential on a supplied list of points ==
      tol = 1d-8
      do i=1,ipol
       call potential_list(ao_bas_han, g_dens(i), npts, qxyz, 
     &                     amat_coul(1,i), tol)
      enddo
      if (ipol.gt.1) then ! enter-ONLY-if ipol>1
       do k = 1,npts
        amat_coul(k,1)=amat_coul(k,1)+amat_coul(k,2)
       enddo
      endif
c
c     == calculate the total nuclear potential on the grid ==
      if (ofinite) then
c ------ Choosing Nuclear Model: erf(zetanuc^0.5 r_L)
        call gridNuclearPotentialFinite(geom,natoms,npts,qxyz,qwght,
     &                             zetanuc_arr,
     &                             closegridpts,amat_nucl)
c ------ Choosing Nuclear Model: P(1/2,zetanuc r_L^2)
c        call gridNuclearPotentialFinite2(geom,natoms,npts,qxyz,qwght,
c     &                             closegridpts,amat_nucl)
      else ! default : point charge model for nuclei
        call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                            closegridpts,amat_nucl)
      endif
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
      if (Knucl) then ! V=Vnucl
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot = amat_nucl(k)
        fac1_arr(k)=amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
     &              *qwght(k)
        denom = (2.d0*clight_au2 - amat_tot)
        fac2_arr(k)=(clight_au2/denom/denom)
     &              *qwght(k)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (do_NonRel) then ! remove it after TEST
          fac1_arr(k)=0.0d0             ! remove it after TEST
          fac2_arr(k)=0.0d0             ! remove it after TEST
        endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
       enddo ! end loop k grid points
      else            ! V=Vnucl+Vee (default)
c     === define fac1_arr,fac2_arr -- FA
       do k = 1,npts
c      == assemble hartree and nuclear contributions ==
        amat_tot =-amat_coul(k,1)+amat_nucl(k)
        fac1_arr(k)=amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
     &              *qwght(k)
        denom = (2.d0*clight_au2 - amat_tot)
        fac2_arr(k)=(clight_au2/denom/denom)
     &              *qwght(k)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (do_NonRel) then ! remove it after TEST
          fac1_arr(k)=0.0d0             ! remove it after TEST
          fac2_arr(k)=0.0d0             ! remove it after TEST
        endif                           ! remove it after TEST
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
       enddo ! end loop k grid points
      endif ! end-if-Knucl
c     == assemble zora correction ==
c -----main diagonal --- START
      do i = 1, nbf
         j=i
         call get_ints_zora_sf(nbf,npts,delchi_ao,i,j,          
     &                         fac1_arr,fac2_arr,      
     &                         ac_sf,ac_scl)  ! out
         zorasf(i,j)   = zorasf(i,j)   + ac_sf
         zorascal(i,j) = zorascal(i,j) + ac_scl
      enddo ! end-loop-i
c -----main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
         call get_ints_zora_sf(nbf,npts,delchi_ao,i,j,          
     &                         fac1_arr,fac2_arr,      
     &                         ac_sf,ac_scl)  ! out
         zorasf(i,j)   = zorasf(i,j)   + 2.0d0*ac_sf
         zorascal(i,j) = zorascal(i,j) + 2.0d0*ac_scl
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
      return
      end
      subroutine get_ints_zora_sf(nbf,            ! in: nr. basis functions
     &                            npts,           ! in: grid points 
     &                            delchi_ao,      ! in: deriv. of basis functions
     &                            i,j,            ! in: (i,j) indices for delchi_ao
     &                            fac1_arr,       ! in
     &                            fac2_arr,       ! in 
     &                            ac_sf,ac_scl)   ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl 
      double precision prod0
          ac_sf    = 0.0d0
          ac_scl   = 0.0d0   
          do k = 1, npts
           prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &            +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &            +delchi_ao(k,3,i)*delchi_ao(k,3,j)  
           ac_sf    = ac_sf    + fac1_arr(k)*prod0
           ac_scl   = ac_scl   + fac2_arr(k)*prod0   
          enddo ! end-loo-k
      return
      end
c $Id$
