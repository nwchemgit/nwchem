c
c     == calculate EFGZ4-scalar relativistic (SR) contribution
      subroutine calc_zora_EFGZ4_SR(ao_bas_han,    ! in: AO basis handle 
     &                              geom,          ! in: geometry handle 
     &                              ipol,          ! in: nr. of polarizations 
     &                              g_dens,        ! in: superposit. atomic densities 
     &                              delchi_ao,     ! in: deriv. of basis functions
     &                              qxyz,          ! in: grid points 
     &                              qwght,         ! in: weighting coeffs.
     &                              nbf,           ! in: nr. basis functions
     &                              npts,          ! in: nr. grid points
     &                              natoms,        ! in: nr. atoms 
     &                              zora_Qpq,      ! in : type of EFG potential
     &                              xyz_EFGcoords, ! in : EFG-nuclear coordinates
     &                              zora_efgz4)    ! out: munu-matrix
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "zora.fh"
#include "global.fh"
#include "bas.fh"

      integer nbf, npts, ao_bas_han, natoms, geom
      integer g_dens(2),ipol
      double precision qwght(npts),pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      double precision amat_Qnucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac_arr(npts)
      integer zora_Qpq
      double precision xyz_EFGcoords(3)      
      double precision ac_efgz4
      double precision zora_efgz4(nbf,nbf)
      external get_ints_zora_efgz4_sr,gridQpqPotential

      clight_au2 = clight_au*clight_au
c
c     == preliminaries ==
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
c     == calculate the total nuclear potential on the grid ==
      call gridNuclearPotentialPoint(geom,natoms,npts,qxyz,qwght,
     &                          closegridpts,amat_nucl)
      call gridQpqPotential(zora_Qpq,xyz_EFGcoords,
     &                      npts,qxyz,           
     &                      amat_Qnucl, ! out: EFG potential
     &                      closegridpts)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac1_arr,fac2_arr -- FA
      do k = 1,npts
c      == assemble hartree and nuclear contributions ==
       amat_tot =-amat_coul(k,1)+amat_nucl(k)
       denom = (2.d0*clight_au2 - amat_tot)
       fac_arr(k)=(clight_au2/denom/denom)*amat_Qnucl(k)
     &             *qwght(k)
      end do
c     == assemble zora correction ==
c -----main diagonal --- START
      do i = 1, nbf
         j=i
         call get_ints_zora_efgz4_sr(nbf,npts,delchi_ao,i,j,          
     &                               fac_arr,  
     &                               ac_efgz4)  ! out
         zora_efgz4(i,j) = zora_efgz4(i,j) + ac_efgz4
      enddo ! end-loop-i
c -----main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
         call get_ints_zora_efgz4_sr(nbf,npts,delchi_ao,i,j,          
     &                               fac_arr,  
     &                               ac_efgz4)  ! out
         zora_efgz4(i,j) = zora_efgz4(i,j)  + 2.0d0*ac_efgz4
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
      return
      end
      subroutine get_ints_zora_efgz4_sr(nbf,            ! in: nr. basis functions
     &                                  npts,           ! in: grid points 
     &                                  delchi_ao,      ! in: deriv. of basis functions
     &                                  i,j,            ! in: (i,j) indices for delchi_ao
     &                                  fac_arr,        ! in
     &                                  ac_efgz4)       ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k
      double precision delchi_ao(npts,3,nbf)
      double precision fac_arr(npts)
      double precision ac_efgz4 
      double precision prod0
          ac_efgz4 = 0.0d0
          do k = 1, npts
           prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &            +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &            +delchi_ao(k,3,i)*delchi_ao(k,3,j)  
           ac_efgz4 = ac_efgz4 + fac_arr(k)*prod0
          enddo ! end-loo-k
      return
      end
c $Id$
