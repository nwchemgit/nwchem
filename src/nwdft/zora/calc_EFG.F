c     == calculate EFG
      subroutine calc_EFG(ao_bas_han,    ! in: AO basis handle 
     &                    geom,          ! in: geometry handle 
     &                    ipol,          ! in: nr. of polarizations 
     &                    g_dens,        ! in: superposit. atomic densities 
     &                    chi_ao,        ! in: basis functions
     &                    qxyz,          ! in: grid points 
     &                    qwght,         ! in: weighting coeffs.
     &                    nbf,           ! in: nr. basis functions
     &                    npts,          ! in: nr. grid points
     &                    natoms,        ! in: nr. atoms 
     &                    zora_Qpq,      ! in : type of EFG potential
     &                    xyz_EFGcoords, ! in : EFG-nuclear coordinates
     &                    zora_efg)    ! out: munu-matrix
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
      double precision chi_ao(npts,nbf)
      integer i,j,k
      double precision amat_Qnucl(npts)
      integer ipt,closegridpts(npts)
      double precision fac_arr(npts)
      integer zora_Qpq
      double precision xyz_EFGcoords(3)      
      double precision ac_efg
      double precision zora_efg(nbf,nbf)
      external get_ints_efg,gridQpqPotential
c     == preliminaries ==
      do ipt = 1,npts
        amat_Qnucl(ipt) = 0.d0
        closegridpts(ipt) = 0
      end do
c     == calculate the total nuclear potential on the grid ==
      call gridQpqPotential(zora_Qpq,xyz_EFGcoords,
     &                      npts,qxyz,           
     &                      amat_Qnucl, ! out: EFG potential
     &                      closegridpts)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac_arr
      do k = 1,npts
       fac_arr(k)=amat_Qnucl(k)*qwght(k)
      end do
c     == assemble zora correction ==
c -----main diagonal --- START
      do i = 1, nbf
         j=i
         call get_ints_efg(nbf,npts,chi_ao,i,j,          
     &                          fac_arr,  
     &                          ac_efg)  ! out
         zora_efg(i,j) = zora_efg(i,j) + ac_efg
      enddo ! end-loop-i
c -----main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
         call get_ints_efg(nbf,npts,chi_ao,i,j,          
     &                          fac_arr,  
     &                          ac_efg)  ! out
         zora_efg(i,j) = zora_efg(i,j)  + 2.0d0*ac_efg
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
      return
      end
      subroutine get_ints_efg(nbf,    ! in: nr. basis functions
     &                        npts,   ! in: grid points 
     &                        chi_ao, ! in: basis functions
     &                        i,j,    ! in: (i,j) indices for delchi_ao
     &                        fac_arr,! in
     &                        ac_efg) ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k
      double precision chi_ao(npts,nbf)
      double precision fac_arr(npts)
      double precision ac_efg
      double precision prod0
          ac_efg = 0.0d0
          do k = 1, npts
           prod0 = chi_ao(k,i)*chi_ao(k,j)
           ac_efg = ac_efg + fac_arr(k)*prod0
          enddo ! end-loo-k
      return
      end
c $Id$
