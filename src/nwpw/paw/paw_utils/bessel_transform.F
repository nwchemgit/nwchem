*
* $Id$
*

!*************************************************
!
!   Name    : spher_bessel_transform
!
!
!   Purpose : performs spherical Bessel transformation
!             from r-space to k-space. The end result is
!             the radial part of the transform in k-space
!
!
!   Created : 8/12/98, modified 7-30-2002
!
!*************************************************
      real*8 function spher_bessel_transform(l,ngrid,log_amesh,r,phi_r,
     >                                       G)
      implicit none
      integer l
      integer ngrid
      real*8 log_amesh
      real*8 r(*)
      real*8 phi_r(*)
      real*8 G

#include "bafdecls.fh"

      !*** Local variables ***
      logical value
      integer i
      integer power_r,power_f
      integer my_sin(2),my_cos(2),tmp_r(2),J(2)
      real*8 phi_G

      !*** external functions ***
#include "integrate.fh"


      value =           BA_push_get(mt_dbl,ngrid,
     >                              'my_sin',my_sin(2),my_sin(1))
      value = value.and.BA_push_get(mt_dbl,ngrid,
     >                              'my_cos',my_cos(2),my_cos(1))
      value = value.and.BA_push_get(mt_dbl,ngrid,
     >                              'tmp_r',tmp_r(2),tmp_r(1))
      value = value.and.BA_push_get(mt_dbl,ngrid,'J',J(2),J(1))
      if (.not.value) 
     >   call errquit('spher_bessel_transform: push stack',0,0)

      do i=1,ngrid
        dbl_mb(tmp_r(1)+i-1)  = r(i)*G
        dbl_mb(my_sin(1)+i-1) = sin(r(i)*G)
        dbl_mb(my_cos(1)+i-1) = cos(r(i)*G)
      end do
   
      power_r = 1
      power_f = 2*l+1

      call find_spher_J_v(l,ngrid,
     >                      dbl_mb(tmp_r(1)), 
     >                      dbl_mb(my_cos(1)), 
     >                      dbl_mb(my_sin(1)),
     >                      dbl_mb(J(1)))

      do i=1,ngrid
        dbl_mb(tmp_r(1)+i-1) = phi_r(i)*dbl_mb(J(1)+i-1)
      end do
      phi_G = def_integr(power_f,dbl_mb(tmp_r(1)),
     >                   power_r,r,
     >                   log_amesh,ngrid)
         
      value =           BA_pop_stack(J(2))
      value = value.and.BA_pop_stack(tmp_r(2))
      value = value.and.BA_pop_stack(my_cos(2))
      value = value.and.BA_pop_stack(my_sin(2)) 
      if (.not.value) call errquit('spher_bessel_transform',0,1)

      spher_bessel_transform = phi_G
      return
      end

!*************************************************
!
!   Name    : spher_bessel0_transform
!
!
!   Purpose : performs spherical Bessel transformation
!             from r-space to k-space. The end result is
!             the radial part of the transform in k-space
!
!
!   Created : 8/12/98, modified 7-30-2002
!
!*************************************************
      real*8 function spher_bessel0_transform(l,ngrid,log_amesh,r,phi_r)
      implicit none
      integer l
      integer ngrid
      real*8 log_amesh
      real*8 r(*)
      real*8 phi_r(*)

      !*** Local variables ***
      real*8 phi_G

      !*** external functions ***
#include "integrate.fh"
      if (l.eq.0) then
        phi_G = def_integr(1,phi_r,
     >                     1,r,
     >                     log_amesh,ngrid)
      else 
        phi_G = 0.0d0
      end if

      spher_bessel0_transform = phi_G
      return
      end

!*************************************************
!
!   Name    : find_spher_J_v
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
      subroutine find_spher_J_v(order,ngrid,x,my_cos,my_sin,J) 
      implicit none
      integer order,ngrid
      real*8  x(*)
      real*8  my_cos(*)
      real*8  my_sin(*)
      real*8  J(*)

*     **** local variables ****
      integer k

      go to (100,200,300,400,500,600,700), order+1
      call errquit('find_spher_J_v: bad order',0,order)

*     *** l=0, s projectors ***
 100  do k=1,ngrid
        J(k) = my_sin(k)/x(k)
      end do
      go to 999

*     *** l=1, p projectors ***
 200  do k=1,ngrid
         J(k)=my_sin(k)/x(k)**2 - my_cos(k)/x(k) 
      end do
      go to 999

*     *** l=2, d projectors ***
 300  do k=1,ngrid
          J(k)=((-3.0d0*my_cos(k))/x(k) 
     >         - my_sin(k) 
     >         + (3.0d0*my_sin(k))/x(k)**2)/x(k)
      end do
      go to 999

*     *** l=3, f projectors ***
 400  do k=1,ngrid
           J(k)=(my_cos(k) 
     >          - (15.0d0*my_cos(k))/x(k)**2 
     >          + (15.0d0*my_sin(k))/x(k)**3 
     >          - ( 6.0d0*my_sin(k))/x(k))/x(k)
      end do
      go to 999

*      *** l=4, g projectors ***
 500  do k=1,ngrid
          J(k)= ((-105.0d0*my_cos(k))/x(k)**3 
     >           + (10.0d0*my_cos(k))/x(k) 
     >           + my_sin(k) 
     >           + (105.0d0*my_sin(k))/x(k)**4
     >           -(45.0d0*my_sin(k))/x(k)**2)/x(k)
      end do
      go to 999

*      *** l=5, h projectors ***
 600  do k=1,ngrid
          J(k)=(-my_cos(k) 
     >         - (945.0d0*my_cos(k))/x(k)**4 
     >         + (105.0d0*my_cos(k))/x(k)**2 
     >         + (945.0d0*my_sin(k))/x(k)**5 
     >         - (420.0d0*my_sin(k))/x(k)**3 
     >         + ( 15.0d0*my_sin(k))/x(k))/x(k)
      end do
      go to 999

*      *** l=6, i projectors ***
 700  do k=1,ngrid
          J(k)=( (-10395.0d0*my_cos(k))/x(k)**5 
     >         + (  1260.0d0*my_cos(k))/x(k)**3 
     >         - (    21.0d0*my_cos(k))/x(k) 
     >         - my_sin(k) 
     >         + (10395.0d0*my_sin(k))/x(k)**6 
     >         - ( 4725.0d0*my_sin(k))/x(k)**4 
     >         + (  210.0d0*my_sin(k))/x(k)**2)/x(k) 
      end do
      go to 999

 999  continue
 
      return
      end !SUBROUTINE find_spher_J_v



