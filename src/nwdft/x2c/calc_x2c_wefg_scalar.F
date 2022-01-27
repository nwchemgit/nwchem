c
c     calculate EFG operator matrix in small component basis, for x2c
      subroutine calc_x2c_wefg_scalar(
     &   ao_bas_han,            ! in: AO basis handle 
     &   geom,                  ! in: geometry handle 
     &   ipol,                  ! in: nr. of polarizations
     &   g_dens,                ! in: superposit. atomic
     &   delchi_ao,             ! in: deriv. basis fcts.
     &   qxyz,                  ! in: grid points 
     &   qwght,                 ! in: weighting coeffs.
     &   nbf,                   ! in: nr. basis functions
     &   npts,                  ! in: nr. grid points
     &   natoms,                ! in: nr. atoms
     &   origin,                ! in: EFG center
     &   x2c0)                  ! out
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
cinclude "x2c.fh"
#include "global.fh"
#include "bas.fh"
c
      integer nbf,npts,ao_bas_han,natoms,geom
      integer g_dens(2),ipol
      double precision qwght(npts)
      double precision pot(npts)
      double precision qxyz(3,npts)
      double precision delchi_ao(npts,3,nbf)
      double precision origin(3)

      double precision x2c0(nbf,nbf,6)

      integer i,j,k, iefg

      integer closegridpts(npts)
      double precision clight_au2
      double precision amat_tot
      double precision denom
      double precision fac1_arr(npts,6)

      double precision ac_sf(6)
  
      external get_ints_x2c_so_efg

c     == preliminaries ==
      do k = 1,npts
        closegridpts(k) = 0
      end do

      do iefg = 1,6
        do k = 1,npts
          fac1_arr(k,iefg) = 0.0d0
        end do
      end do

c     == calculate the EFG operator on the grid ==
            
        call gridQpqPotential_all(origin,
     &     npts,qxyz,
     &     fac1_arr,            ! out: EFG potential
     &     closegridpts)

c     black out grid points too close to a nucleus;
c     the numerical integration may get noisy otherwise:

      do k = 1,npts
        if (k.eq.closegridpts(k)) qwght(k) = 0.d0
      end do

c     === premultiply fac1_arr with integration weights
      do iefg = 1,6
        do k = 1,npts
          fac1_arr(k,iefg)= fac1_arr(k,iefg) * qwght(k)
        end do
      end do

c     == assemble integrals ==

c     ----- main diagonal --- START
      do i = 1, nbf
        j=i
        call get_ints_x2c_sc_efg(
     &     nbf,npts,delchi_ao,i,j,      
     &     fac1_arr,
     &     ac_sf)             ! out 

        do iefg = 1,6
          x2c0(i,j,iefg) = x2c0(i,j,iefg) + ac_sf(iefg)
        end do
      enddo                     ! end-loop-i
c     ----- main diagonal --- END
      
c     ----- off diagonal --- START
      do i = 1, nbf
        do j = i+1, nbf
          call get_ints_x2c_sc_efg(
     &       nbf,npts,delchi_ao,i,j,      
     &       fac1_arr,
     &       ac_sf)             ! out 

          do iefg = 1,6
            x2c0(i,j,iefg) = x2c0(i,j,iefg) + 2.0d0*ac_sf(iefg)
          end do
        enddo                   ! end-loop-j
      enddo                     ! end-loop-i
c     ----- off diagonal --- END

c     all done
      return
      end


      subroutine get_ints_x2c_sc_efg(
     &   nbf,                   ! in: # basis fncts.
     &   npts,                  ! in: grid points 
     &   delchi_ao,             ! in: deriv. of basis fcts.
     &   i,j,                   ! in: (i,j) indices for delchi_ao
     &   fac1_arr,              ! operator on grid points * weights
     &   ac_sf)                 ! out


      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer nbf,npts,i,j,k, iefg
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts,6)
      double precision
     &   ac_sf(6)
      double precision prod0

      do iefg = 1,6
        ac_sf(iefg)    = 0.0d0
      end do
      
      do k = 1, npts

        prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &         +delchi_ao(k,2,i)*delchi_ao(k,2,j) 
     &         +delchi_ao(k,3,i)*delchi_ao(k,3,j)  

        do iefg = 1,6
          ac_sf(iefg)    = ac_sf(iefg) + fac1_arr(k,iefg)*prod0
        end do                  ! iefg
        
      end do                    ! k

      return
      end

