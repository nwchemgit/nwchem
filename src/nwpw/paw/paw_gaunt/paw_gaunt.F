c
c $Id: paw_gaunt.F,v 1.5 2004-07-29 15:55:37 bylaska Exp $
c
c $Log: not supported by cvs2svn $
c Revision 1.4  2003/10/21 02:05:16  marat
c switched to new errquit by running global replace operation
c see the script below (note it will not work on multiline errquit calls)
c *********************************************************
c #!/bin/sh
c
c e=`find . -name "*F" -print`
c
c for f in $e
c do
c cp $f $f.bak
c sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
c #rm $f.bak
c done
c **********************************************************
c
c Revision 1.3  2003/02/11 01:41:23  edo
c eliminated f90-isms
c
c Revision 1.2  2003/02/10 02:36:26  marat
c ...
c
c Revision 1.1  2003/02/06 06:04:37  marat
c Initial Revision
c

*     ***********************************
*     *                                 *
*     *  	 paw_gaunt_init 	*
*     *                                 *
*     ***********************************
      subroutine paw_gaunt_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"

#include "paw_gaunt_data.fh"

*     !**** local varialbles ****
      logical value
      integer i,j,k,indx
      integer l1,l2,l3,m1,m2,m3

      integer sizel,sizel2

      !**** external functions ****
      real*8   gen_gaunt_coeff_sub
      external gen_gaunt_coeff_sub

      call nwpw_timing_start(15)

      paw_gaunt_lmax = 8

      sizel  = (paw_gaunt_lmax+1)
      sizel2 = sizel**2

      value = MA_alloc_get(mt_dbl,sizel*sizel2*sizel2,
     >                    'paw_gaunt_coeff',
     >                     paw_gaunt_coeff(2),
     >                     paw_gaunt_coeff(1))
      if (.not. value) 
     >   call errquit('paw_gaunt_init:out of heap memory',0,
     >                MA_ERR)


      do l2=0,paw_gaunt_lmax
      do l3=0,paw_gaunt_lmax

        do m2=-l2,l2
        do m3=-l3,l3
          m1 = m2-m3

          do l1=0,abs(m1)
            i = l1
            j = l2**2 + (l2+m2)
            k = l3**2 + (l3+m3)
            indx = i + j*sizel + k*sizel*sizel2

            dbl_mb(paw_gaunt_coeff(1)+indx)
     >      =gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3) 
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
*     *  	 paw_gaunt_end  	*
*     *                                 *
*     ***********************************
      subroutine paw_gaunt_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"

#include "paw_gaunt_data.fh"

      logical value

      call nwpw_timing_start(15)

      value = MA_free_heap(paw_gaunt_coeff(2)) 
      if (.not. value) 
     > call errquit('paw_gaunt_end:error freeing heap memory',0,MA_ERR)

      call nwpw_timing_end(15)
      return
      end

*     ***********************************
*     *                                 *
*     *  	 gen_gaunt_coeff	*
*     *                                 *
*     ***********************************
      double precision function gen_gaunt_coeff(l1,m1,l2,m2,l3,m3) 
      implicit none
      integer  l1,l2,l3
      integer  m1,m2,m3

#include "mafdecls.fh"
#include "errquit.fh"

#include "paw_gaunt_data.fh"


      !**** local variables ****
      integer i,j,k,indx
      integer sizel,sizel2
      real*8 sum

      !**** external functions ****
      real*8   gen_gaunt_coeff_sub
      external gen_gaunt_coeff_sub

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
      else if ((l1.gt.paw_gaunt_lmax) .or. 
     >         (l2.gt.paw_gaunt_lmax) .or. 
     >         (l3.gt.paw_gaunt_lmax)) then
        sum = gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3)

      !*** retrieve precomputed gaunt coefficient ***
      else
        sizel  = (paw_gaunt_lmax+1)
        sizel2 = sizel**2
        i = l1
        j = l2**2 + (l2+m2)
        k = l3**2 + (l3+m3)
        indx = i + j*sizel + k*sizel*sizel2
        sum = dbl_mb(paw_gaunt_coeff(1)+indx)
      end if

      call nwpw_timing_end(15)

      gen_gaunt_coeff = sum
      return
      end


*     ***********************************
*     *                                 *
*     *      gen_gaunt_coeff_sub	*
*     *                                 *
*     ***********************************

      real*8 function gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3)
      implicit none
      integer l1, l2, l3
      integer m1, m2, m3
      integer i
      integer order
      double precision x1, x2
     

#include"paw_params.fh"
#include"paw_spher_func.fh"

c     work arrays for integration
      double precision x(100),w(100)


      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) call errquit(
     1   'Invalid parameter in gen_gaunt_coeff, negative l',0,1)
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2)) call 
     1   errquit('Invalid parameter in gen_gaunt_coeff, m > l',0,1)

      !**** Do integration over angle phi ****
      if ((-m1) + m2 - m3 .ne. 0) then
         gen_gaunt_coeff_sub = 0.0d0
         return 
      endif

      !**** Check the triangle rule ****
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         gen_gaunt_coeff_sub = 0.0d0
         return 
      endif

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         gen_gaunt_coeff_sub = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         gen_gaunt_coeff_sub = 1.0d0/dsqrt(fourpi)
         return 
      endif
      x1 = -1.0
      x2 = 1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL gauss_weights(x1,x2,x,w,order)
      gen_gaunt_coeff_sub = 0.0d0
      do i = 1, order
         gen_gaunt_coeff_sub = gen_gaunt_coeff_sub 
     >                       + w(i)*theta_lm(l1,m1,x(i))*
     >      theta_lm(l2,m2,x(i))*theta_lm(l3,m3,x(i))
      end do

      gen_gaunt_coeff_sub = gen_gaunt_coeff_sub*2.0d0*pi

      return
      end 
 
