c
c $Id$
c

*     ***********************************
*     *                                 *
*     *  	 nwpw_gaunt_init 	*
*     *                                 *
*     ***********************************
      subroutine paw_gaunt_init(iscomplex)
      implicit none
      logical iscomplex

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_gaunt_data.fh"

*     !**** local varialbles ****
      logical value
      integer i,j,k,indx
      integer l1,l2,l3,m1,m2,m3

      integer sizel,sizel2

      !**** external functions ****
      real*8   nwpw_gen_gaunt_coeff_sub
      external nwpw_gen_gaunt_coeff_sub

      call nwpw_timing_start(15)

      gaunt_complex = iscomplex
      gaunt_lmax = 8

      sizel  = (gaunt_lmax+1)
      sizel2 = sizel**2

      value = BA_alloc_get(mt_dbl,sizel*sizel2*sizel2,
     >                    'gaunt_coeff',
     >                     gaunt_coeff(2),
     >                     gaunt_coeff(1))
      if (.not. value) 
     >   call errquit('nwpw_gaunt_init:out of heap memory',0,
     >                MA_ERR)


      do l1=0,gaunt_lmax
      do l2=0,gaunt_lmax
      do l3=0,gaunt_lmax

        do m2=-l2,l2
        do m3=-l3,l3
          m1 = m2-m3

          if (abs(m1) .le. l1) then
            i = l1
            j = l2**2 + (l2+m2)
            k = l3**2 + (l3+m3)
            indx = i + j*sizel + k*sizel*sizel2

            dbl_mb(gaunt_coeff(1)+indx)
     >      =nwpw_gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3) 
          end if

        end do
        end do

      end do
      end do
      end do
      
      call nwpw_timing_end(15)
      return
      end

*     ***********************************
*     *                                 *
*     *  	 nwpw_gaunt_end  	*
*     *                                 *
*     ***********************************
      subroutine nwpw_gaunt_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_gaunt_data.fh"


      call nwpw_timing_start(15)

      if (.not.BA_free_heap(nwpw_gaunt_coeff(2)))
     > call errquit('paw_gaunt_end:error freeing heap memory',0,MA_ERR)

      call nwpw_timing_end(15)
      return
      end

*     ***********************************
*     *                                 *
*     *       nwpw_gen_gaunt_coeff	*
*     *                                 *
*     ***********************************
      double precision function nwpw_gen_gaunt_coeff(l1,m1,l2,m2,l3,m3)
      implicit none
      integer  l1,l2,l3
      integer  m1,m2,m3

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_gaunt_data.fh"


      !**** local variables ****
      integer i,j,k,indx
      integer sizel,sizel2
      real*8 sum

      !**** external functions ****
      real*8   nwpw_gen_gaunt_coeff_sub
      external nwpw_gen_gaunt_coeff_sub

      call nwpw_timing_start(15)

!      IF(l1<0 .OR. l2<0 .OR. l3<0) THEN
!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, negative l")
!      END IF
!
!      IF(l1<ABS(m1) .OR. l3<ABS(m3) .OR. l2<ABS(m2)) THEN
!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, m > l")
!      END IF

      !*** gaunt coefficient zero ***
      if ((-m1+m2-m3).ne.0) then
        sum = 0.0d0

      !*** calculate gaunt coefficient ***
      else if ((l1.gt.gaunt_lmax) .or. 
     >         (l2.gt.gaunt_lmax) .or. 
     >         (l3.gt.gaunt_lmax)) then
        sum = nwpw_gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3)

      !*** retrieve precomputed gaunt coefficient ***
      else
        sizel  = (gaunt_lmax+1)
        sizel2 = sizel**2
        i = l1
        j = l2**2 + (l2+m2)
        k = l3**2 + (l3+m3)
        indx = i + j*sizel + k*sizel*sizel2
        sum = dbl_mb(gaunt_coeff(1)+indx)
      end if

      call nwpw_timing_end(15)

      nwpw_gen_gaunt_coeff = sum
      return
      end


*     ***********************************
*     *                                 *
*     *      nwpw_gen_gaunt_coeff_sub	*
*     *                                 *
*     ***********************************

      real*8 function nwpw_gen_gaunt_coeff_sub(iscmplx,
     >                                         l1,m1,l2,m2,l3,m3)
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     

#include"paw_spher_func.fh"


*     *** local variables ****
      double precision pi
      double precision fourpi
      parameter(pi = 3.14159265358979323846264338327950288419D0)
      parameter(fourpi = pi*4.0D0)

c     work arrays for integration
      double precision x(100),w(100)


      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) call errquit(
     1   'Invalid parameter in gen_gaunt_coeff, negative l',0,1)
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2)) call 
     1   errquit('Invalid parameter in gen_gaunt_coeff, m > l',0,1)

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gen_gaunt_coeff_sub = 0.0d0
            return 
         endif
      else
      endif

      !**** Check the triangle rule ****
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         nwpw_gen_gaunt_coeff_sub = 0.0d0
         return 
      endif

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         nwpw_gen_gaunt_coeff_sub = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gen_gaunt_coeff_sub = 1.0d0/dsqrt(fourpi)
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL gauss_weights(x1,x2,x,w,order)
      nwpw_gen_gaunt_coeff_sub = 0.0d0
      do i = 1, order
         nwpw_gen_gaunt_coeff_sub = nwpw_gen_gaunt_coeff_sub 
     >                       + w(i)*theta_lm(l1,m1,x(i))
     >                             *theta_lm(l2,m2,x(i))
     >                             *theta_lm(l3,m3,x(i))
      end do

      nwpw_gen_gaunt_coeff_sub = nwpw_gen_gaunt_coeff_sub*2.0d0*pi

      return
      end 
 
