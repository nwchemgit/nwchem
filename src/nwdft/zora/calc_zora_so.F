c
c     == calculate spin-orbit zora contribution ==
      subroutine calc_zora_so(ao_bas_han,  ! in: AO basis handle 
     &                           geom,        ! in: geometry handle 
     &                           ipol,        ! in: nr. of polarizations 
     &                           g_dens,      ! in: superposit. atomic densities 
     &                           delchi_ao,   ! in: deriv. of basis functions
     &                           qxyz,        ! in: grid points 
     &                           qwght,       ! in: weighting coeffs.
     &                           nbf,         ! in: nr. basis functions
     &                           npts,        ! in: nr. grid points
     &                           natoms,      ! in: nr. atoms 
     &                           zora0, scal0,! out  
     &                           zorax, scalx,! out 
     &                           zoray, scaly,! out 
     &                           zoraz, scalz)! out
      implicit none
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
      double precision delchi_ao(npts,3,nbf)
      double precision zora0(nbf,nbf)  ! scalar component
      double precision scal0(nbf,nbf)  ! scalar component
      double precision zorax(nbf,nbf)
      double precision scalx(nbf,nbf)
      double precision zoray(nbf,nbf)
      double precision scaly(nbf,nbf)
      double precision zoraz(nbf,nbf)
      double precision scalz(nbf,nbf)
      integer i,j,k
      double precision amat_coul(npts,ipol)
      double precision amat_nucl(npts)
      integer ipt,closegridpts(npts)
      double precision clight_au2,tol
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl,
     &                 ac_sx,ac_scl_x,
     &                 ac_sy,ac_scl_y,
     &                 ac_sz,ac_scl_z      
      external get_ints_zora_so

      clight_au2 = clight_au*clight_au
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
      call gridNuclearPotential(geom,natoms,npts,qxyz,qwght,
     &                          closegridpts,amat_nucl)
      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do
c     === define fac1_arr,fac2_arr -- FA
      do k = 1,npts
c      == assemble hartree and nuclear contributions ==
       amat_tot =-amat_coul(k,1)+amat_nucl(k)
       fac1_arr(k)=amat_tot/(4.d0*clight_au2 - 2.d0*amat_tot)
     &             *qwght(k)
       denom = (2.d0*clight_au2 - amat_tot)
       fac2_arr(k)=(clight_au2/denom/denom)
     &             *qwght(k)
      end do
c     == assemble zora correction ==
c ----- main diagonal --- START
      do i = 1, nbf
          j=i
          call get_ints_zora_so(nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,fac2_arr,       
     &                          ac_sf,ac_scl,   ! out
     &                          ac_sx,ac_scl_x, ! out 
     &                          ac_sy,ac_scl_y, ! out 
     &                          ac_sz,ac_scl_z) ! out 
          zora0(i,j) = zora0(i,j) + ac_sf
          scal0(i,j) = scal0(i,j) + ac_scl
          zorax(i,j) = zorax(i,j) + ac_sx
          scalx(i,j) = scalx(i,j) + ac_scl_x
          zoray(i,j) = zoray(i,j) + ac_sy
          scaly(i,j) = scaly(i,j) + ac_scl_y
          zoraz(i,j) = zoraz(i,j) + ac_sz
          scalz(i,j) = scalz(i,j) + ac_scl_z
      enddo ! end-loop-i
c ----- main diagonal --- END
c ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_zora_so(nbf,npts,delchi_ao,i,j,      
     &                          fac1_arr,fac2_arr,       
     &                          ac_sf,ac_scl,   ! out
     &                          ac_sx,ac_scl_x, ! out 
     &                          ac_sy,ac_scl_y, ! out 
     &                          ac_sz,ac_scl_z) ! out 
          zora0(i,j) = zora0(i,j) + 2.0d0*ac_sf
          scal0(i,j) = scal0(i,j) + 2.0d0*ac_scl
          zorax(i,j) = zorax(i,j) + 2.0d0*ac_sx
          scalx(i,j) = scalx(i,j) + 2.0d0*ac_scl_x
          zoray(i,j) = zoray(i,j) + 2.0d0*ac_sy
          scaly(i,j) = scaly(i,j) + 2.0d0*ac_scl_y
          zoraz(i,j) = zoraz(i,j) + 2.0d0*ac_sz
          scalz(i,j) = scalz(i,j) + 2.0d0*ac_scl_z
        enddo ! end-loop-j
      enddo ! end-loop-i
c ----- off diagonal --- END
      return
      end
      subroutine get_ints_zora_so(nbf,            ! in: nr. basis functions
     &                            npts,           ! in: grid points 
     &                            delchi_ao,      ! in: deriv. of basis functions
     &                            i,j,            ! in: (i,j) indices for delchi_ao
     &                            fac1_arr,       ! in
     &                            fac2_arr,       ! in 
     &                            ac_sf,ac_scl,   ! out
     &                            ac_sx,ac_scl_x, ! out 
     &                            ac_sy,ac_scl_y, ! out 
     &                            ac_sz,ac_scl_z) ! out
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl,
     &                 ac_sx,ac_scl_x,
     &                 ac_sy,ac_scl_y,
     &                 ac_sz,ac_scl_z      
      double precision prod0,prod_x,prod_y,prod_z
          ac_sf    = 0.0d0
          ac_sx    = 0.0d0
          ac_sy    = 0.0d0
          ac_sz    = 0.0d0
          ac_scl   = 0.0d0   
          ac_scl_x = 0.0d0
          ac_scl_y = 0.0d0
          ac_scl_z = 0.0d0       
          do k = 1, npts
           prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &            +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &            +delchi_ao(k,3,i)*delchi_ao(k,3,j)  
           prod_x= delchi_ao(k,2,i)*delchi_ao(k,3,j)
     &            -delchi_ao(k,3,i)*delchi_ao(k,2,j) 
           prod_y= delchi_ao(k,3,i)*delchi_ao(k,1,j)
     &            -delchi_ao(k,1,i)*delchi_ao(k,3,j) 
           prod_z= delchi_ao(k,1,i)*delchi_ao(k,2,j)
     &            -delchi_ao(k,2,i)*delchi_ao(k,1,j) 
           ac_sf    = ac_sf    + fac1_arr(k)*prod0
           ac_scl   = ac_scl   + fac2_arr(k)*prod0   
           ac_sx    = ac_sx    + fac1_arr(k)*prod_x
           ac_scl_x = ac_scl_x + fac2_arr(k)*prod_x  
           ac_sy    = ac_sy    + fac1_arr(k)*prod_y
           ac_scl_y = ac_scl_y + fac2_arr(k)*prod_y    
           ac_sz    = ac_sz    + fac1_arr(k)*prod_z
           ac_scl_z = ac_scl_z + fac2_arr(k)*prod_z        
          enddo ! end-loo-k
      return
      end
c $Id$
