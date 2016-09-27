c
      subroutine get_ints_zora(do_so, 
     &                         nbf,             ! in: nr. basis functions
     &                         npts,            ! in: grid points 
     &                         delchi_ao,       ! in: deriv. of basis functions
     &                         i,j,             ! in: (i,j) indices for delchi_ao
     &                         fac1_arr,        ! in
     &                         fac2_arr,        ! in 
     &                         ac_sf,ac_scl,    ! out
     &                         ac_sx,ac_scl_x,  ! out 
     &                         ac_sy,ac_scl_y,  ! out 
     &                         ac_sz,ac_scl_z)  ! out
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
c
      logical do_so
      integer nbf,npts,i,j
      double precision delchi_ao(npts,3,nbf)
      double precision fac1_arr(npts),fac2_arr(npts)
      double precision ac_sf,ac_scl,
     &                 ac_sx,ac_scl_x,
     &                 ac_sy,ac_scl_y,
     &                 ac_sz,ac_scl_z      
c
      integer k
      double precision prod0,prod_x,prod_y,prod_z
c
      ac_sf    = 0.0d0
      ac_scl   = 0.0d0   
      ac_sx    = 0.0d0
      ac_scl_x = 0.0d0
      ac_sy    = 0.0d0
      ac_scl_y = 0.0d0
      ac_sz    = 0.0d0
      ac_scl_z = 0.0d0       
c
      if (do_so) then  ! for spin-orbit
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
        end do ! npts
      else  
        do k = 1, npts
           prod0 = delchi_ao(k,1,i)*delchi_ao(k,1,j)
     &            +delchi_ao(k,2,i)*delchi_ao(k,2,j)
     &            +delchi_ao(k,3,i)*delchi_ao(k,3,j)
           ac_sf    = ac_sf    + fac1_arr(k)*prod0
           ac_scl   = ac_scl   + fac2_arr(k)*prod0
        end do ! npts
      end if  ! do_so
c
      return
      end
