c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.6  2004/07/31 21:54:30  bylaska
c Added templates for paw_energy and paw_gradient routines ...EJB
c
c Revision 1.5  2004/07/29 15:55:37  bylaska
c Temporary array for Gaunt coefficients added.  Speeds up the program considerably, but it is extrememely memory intensive.  Also added timing routines to multipole calculations and fixed the initial total charge calculation.
c
c  ...EJB
c
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
*     *          paw_gaunt_init         *
*     *                                 *
*     ***********************************
      subroutine paw_gaunt_init()
      implicit none

#include "bafdecls.fh"
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

      value = BA_alloc_get(mt_dbl,sizel*sizel2*sizel2,
     >                    'paw_gaunt_coeff',
     >                     paw_gaunt_coeff(2),
     >                     paw_gaunt_coeff(1))
      if (.not. value) 
     >   call errquit('paw_gaunt_init:out of heap memory',0,
     >                MA_ERR)


      do l1=0,paw_gaunt_lmax
      do l2=0,paw_gaunt_lmax
      do l3=0,paw_gaunt_lmax

        do m2=-l2,l2
        do m3=-l3,l3
          m1 = m2-m3

          if (abs(m1) .le. l1) then
            i = l1
            j = l2**2 + (l2+m2)
            k = l3**2 + (l3+m3)
            indx = i + j*sizel + k*sizel*sizel2

            dbl_mb(paw_gaunt_coeff(1)+indx)
     >      =gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3) 
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
*     *          paw_gaunt_end          *
*     *                                 *
*     ***********************************
      subroutine paw_gaunt_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "paw_gaunt_data.fh"

      logical value

      call nwpw_timing_start(15)

      value = BA_free_heap(paw_gaunt_coeff(2)) 
      if (.not. value) 
     > call errquit('paw_gaunt_end:error freeing heap memory',0,MA_ERR)

      call nwpw_timing_end(15)
      return
      end

*     ***********************************
*     *                                 *
*     *          gen_gaunt_coeff        *
*     *                                 *
*     ***********************************
      double precision function gen_gaunt_coeff(l1,m1,l2,m2,l3,m3) 
      implicit none
      integer  l1,l2,l3
      integer  m1,m2,m3

#include "bafdecls.fh"
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
*     *      gen_gaunt_coeff_sub        *
*     *                                 *
*     ***********************************

      real*8 function gen_gaunt_coeff_sub(l1,m1,l2,m2,l3,m3)
      implicit none
      integer l1, l2, l3
      integer m1, m2, m3
      integer i
      integer order
      double precision x1, x2
     

#include "paw_params.fh"
#include "paw_spher_func.fh"

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
 


c
c $Id$
c
c $Log: not supported by cvs2svn $

C*     ***********************************
C*     *                                 *
C*     *         paw_gaunt2_init        *
C*     *                                 *
C*     ***********************************
C*
C*  needed for d/dtheta derivatives in paw_xc
C
C      subroutine paw_gaunt2_init()
C      implicit none
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C*     !**** local varialbles ****
C      logical value
C      integer i,j,k,indx
C      integer l1,l2,l3,m1,m2,m3
C
C      integer sizel,sizel2
C
C      !**** external functions ****
C      real*8   gen_gaunt2_coeff_sub
C      external gen_gaunt2_coeff_sub
C
C      call nwpw_timing_start(15)
C
C      paw_gaunt2_lmax = 8
C
C      sizel  = (paw_gaunt2_lmax+1)
C      sizel2 = sizel**2
C
C*     **** such a large array is probably not needed ****
C      value = BA_alloc_get(mt_dbl,sizel*sizel2*sizel2,
C     >                    'paw_gaunt2_coeff',
C     >                     paw_gaunt2_coeff(2),
C     >                     paw_gaunt2_coeff(1))
C      if (.not. value) 
C     >   call errquit('paw_gaunt_init:out of heap memory',0,
C     >                MA_ERR)
C
C
C      do l1=0,paw_gaunt2_lmax
C      do l2=0,paw_gaunt2_lmax
C      do l3=0,paw_gaunt2_lmax
C
C        do m2=-l2,l2
C        do m3=-l3,l3
C          m1 = m2-m3
C
C          if (abs(m1) .le. l1) then
C            i = l1
C            j = l2**2 + (l2+m2)
C            k = l3**2 + (l3+m3)
C            indx = i + j*sizel + k*sizel*sizel2
C
C            dbl_mb(paw_gaunt2_coeff(1)+indx)
C     >      =gen_gaunt2_coeff_sub(l1,m1,l2,m2,l3,m3) 
C          end if
C
C        end do
C        end do
C
C      end do
C      end do
C      end do
C      
C      call nwpw_timing_end(15)
C      return
C      end
C
C*     ***********************************
C*     *                                 *
C*     *         paw_gaunt2_end         *
C*     *                                 *
C*     ***********************************
C      subroutine paw_gaunt2_end()
C      implicit none
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C      logical value
C
C      call nwpw_timing_start(15)
C
C      value = BA_free_heap(paw_gaunt2_coeff(2)) 
C      if (.not. value) 
C     > call errquit('paw_gaunt_end:error freeing heap memory',0,MA_ERR)
C
C      call nwpw_timing_end(15)
C      return
C      end
C
C*     ***********************************
C*     *                                 *
C*     *         gen_gaunt2_coeff       *
C*     *                                 *
C*     ***********************************
C      double precision function gen_gaunt2_coeff(l1,m1,l2,m2,l3,m3) 
C      implicit none
C      integer  l1,l2,l3
C      integer  m1,m2,m3
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C
C      !**** local variables ****
C      integer i,j,k,indx
C      integer sizel,sizel2
C      real*8 sum
C
C      !**** external functions ****
C      real*8   gen_gaunt2_coeff_sub
C      external gen_gaunt2_coeff_sub
C
C      call nwpw_timing_start(15)
C
C!      IF(l1<0 .OR. l2<0 .OR. l3<0) THEN
C!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, negative l")
C!      END IF
C!
C!      IF(l1<ABS(m1) .OR. l3<ABS(m3) .OR. l2<ABS(m2)) THEN
C!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, m > l")
C!      END IF
C
C      !*** gaunt coefficient zero ***
C      if ((-m1+m2-m3).ne.0) then
C        sum = 0.0d0
C
C      !*** calculate gaunt coefficient ***
C      else if ((l1.gt.paw_gaunt2_lmax) .or. 
C     >         (l2.gt.paw_gaunt2_lmax) .or. 
C     >         (l3.gt.paw_gaunt2_lmax)) then
C        sum = gen_gaunt2_coeff_sub(l1,m1,l2,m2,l3,m3)
C
C      !*** retrieve precomputed gaunt coefficient ***
C      else
C        sizel  = (paw_gaunt2_lmax+1)
C        sizel2 = sizel**2
C        i = l1
C        j = l2**2 + (l2+m2)
C        k = l3**2 + (l3+m3)
C        indx = i + j*sizel + k*sizel*sizel2
C        sum = dbl_mb(paw_gaunt2_coeff(1)+indx)
C      end if
C
C      call nwpw_timing_end(15)
C
C      gen_gaunt2_coeff = sum
C      return
C      end


*     ***********************************
*     *                                 *
*     *      gen_gaunt2_coeff_sub       *
*     *                                 *
*     ***********************************

      real*8 function gen_gaunt2_coeff_sub(l1,m1,l2,m2,l3,m3)
      implicit none
      integer l1, l2, l3
      integer m1, m2, m3
      integer i
      integer order
     

#include "paw_params.fh"
#include "paw_spher_func.fh"

c     work arrays for integration
      double precision x(100),w(100)


      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) call errquit(
     1   'Invalid parameter in gen_gaunt_coeff, negative l',0,1)
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2)) call 
     1   errquit('Invalid parameter in gen_gaunt_coeff, m > l',0,1)

      !**** integration over angle phi = 0 or 2*pi ****
      if ((-m1) + m2 - m3 .ne. 0) then
         gen_gaunt2_coeff_sub = 0.0d0
         return 
      endif

      !**** no triangle rule ****

      !**** Check if the integrand is odd function==>integral is zero ****
      !*** derivative wrt theta makes an reverse odd/even behavior    ****
      if (mod(l1 + l2 + l3,2) .eq. 0) then
         gen_gaunt2_coeff_sub = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l2.eq.0 .and. l3.eq.0) then
         gen_gaunt2_coeff_sub = 0.0d0
         return 
      endif
      order = (l1 + l2 + l3) + 1

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL gauss_chebyshev_weights(x,w,order)
      gen_gaunt2_coeff_sub = 0.0d0
      do i = 1, order
         gen_gaunt2_coeff_sub 
     >    = gen_gaunt2_coeff_sub 
     >    + w(i)*theta_lm(l1,m1,x(i))
     >          *dsqrt((1.0d0-x(i))*(1.0d0+x(i)))
     >          *(dtheta_lm(l2,m2,x(i))* theta_lm(l3,m3,x(i))
     >           + theta_lm(l2,m2,x(i))*dtheta_lm(l3,m3,x(i)))
      end do

      gen_gaunt2_coeff_sub = gen_gaunt2_coeff_sub*2.0d0*pi

      return
      end 
 



C*     ***********************************
C*     *                                 *
C*     *         paw_gaunt3_init        *
C*     *                                 *
C*     ***********************************
C*
C*  needed for 1/sin(theta)*d/dphi derivatives in paw_xc
C
C      subroutine paw_gaunt3_init()
C      implicit none
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C*     !**** local varialbles ****
C      logical value
C      integer i,j,k,indx
C      integer l1,l2,l3,m1,m2,m3
C
C      integer sizel,sizel2
C
C      !**** external functions ****
C      real*8   gen_gaunt3_coeff_sub
C      external gen_gaunt3_coeff_sub
C
C      call nwpw_timing_start(15)
C
C      paw_gaunt3_lmax = 8
C
C      sizel  = (paw_gaunt3_lmax+1)
C      sizel2 = sizel**2
C
C*     **** such a large array is probably not needed ****
C      value = BA_alloc_get(mt_dbl,sizel*sizel2*sizel2,
C     >                    'paw_gaunt3_coeff',
C     >                     paw_gaunt3_coeff(2),
C     >                     paw_gaunt3_coeff(1))
C      if (.not. value) 
C     >   call errquit('paw_gaunt_init:out of heap memory',0,
C     >                MA_ERR)
C
C
C      do l1=0,paw_gaunt3_lmax
C      do l2=0,paw_gaunt3_lmax
C      do l3=0,paw_gaunt3_lmax
C
C        do m2=-l2,l2
C        do m3=-l3,l3
C          m1 = m2-m3
C
C          if (abs(m1) .le. l1) then
C            i = l1
C            j = l2**2 + (l2+m2)
C            k = l3**2 + (l3+m3)
C            indx = i + j*sizel + k*sizel*sizel2
C
C            dbl_mb(paw_gaunt3_coeff(1)+indx)
C     >      =gen_gaunt3_coeff_sub(l1,m1,l2,m2,l3,m3) 
C          end if
C
C        end do
C        end do
C
C      end do
C      end do
C      end do
C      
C      call nwpw_timing_end(15)
C      return
C      end
C
C*     ***********************************
C*     *                                 *
C*     *         paw_gaunt3_end         *
C*     *                                 *
C*     ***********************************
C      subroutine paw_gaunt3_end()
C      implicit none
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C      logical value
C
C      call nwpw_timing_start(15)
C
C      value = BA_free_heap(paw_gaunt3_coeff(2)) 
C      if (.not. value) 
C     > call errquit('paw_gaunt_end:error freeing heap memory',0,MA_ERR)
C
C      call nwpw_timing_end(15)
C      return
C      end
C
C*     ***********************************
C*     *                                 *
C*     *         gen_gaunt3_coeff       *
C*     *                                 *
C*     ***********************************
C      double precision function gen_gaunt3_coeff(l1,m1,l2,m2,l3,m3) 
C      implicit none
C      integer  l1,l2,l3
C      integer  m1,m2,m3
C
C#include "bafdecls.fh"
C#include "errquit.fh"
C
C#include "paw_gaunt_data.fh"
C
C
C      !**** local variables ****
C      integer i,j,k,indx
C      integer sizel,sizel2
C      real*8 sum
C
C      !**** external functions ****
C      real*8   gen_gaunt3_coeff_sub
C      external gen_gaunt3_coeff_sub
C
C      call nwpw_timing_start(15)
C
C!      IF(l1<0 .OR. l2<0 .OR. l3<0) THEN
C!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, negative l")
C!      END IF
C!
C!      IF(l1<ABS(m1) .OR. l3<ABS(m3) .OR. l2<ABS(m2)) THEN
C!       CALL fatal_error("Invalid parameter in gen_gaunt_coeff, m > l")
C!      END IF
C
C      !*** gaunt coefficient zero ***
C      if ((-m1+m2-m3).ne.0) then
C        sum = 0.0d0
C
C      !*** calculate gaunt coefficient ***
C      else if ((l1.gt.paw_gaunt3_lmax) .or. 
C     >         (l2.gt.paw_gaunt3_lmax) .or. 
C     >         (l3.gt.paw_gaunt3_lmax)) then
C        sum = gen_gaunt3_coeff_sub(l1,m1,l2,m2,l3,m3)
C
C      !*** retrieve precomputed gaunt coefficient ***
C      else
C        sizel  = (paw_gaunt3_lmax+1)
C        sizel2 = sizel**2
C        i = l1
C        j = l2**2 + (l2+m2)
C        k = l3**2 + (l3+m3)
C        indx = i + j*sizel + k*sizel*sizel2
C        sum = dbl_mb(paw_gaunt3_coeff(1)+indx)
C      end if
C
C      call nwpw_timing_end(15)
C
C      gen_gaunt3_coeff = sum
C      return
C      end


*     ***********************************
*     *                                 *
*     *      gen_gaunt3_coeff_sub       *
*     *                                 *
*     ***********************************

      real*8 function gen_gaunt3_coeff_sub(l1,m1,l2,m2,l3,m3)
      implicit none
      integer l1, l2, l3
      integer m1, m2, m3
      integer i
      integer order
     

#include "paw_params.fh"
#include "paw_spher_func.fh"

c     work arrays for integration
      double precision x(100),w(100),f


      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) call errquit(
     1   'Invalid parameter in gen_gaunt_coeff, negative l',0,1)
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2)) call 
     1   errquit('Invalid parameter in gen_gaunt_coeff, m > l',0,1)

      !**** integration over angle phi = 0 or 2*pi ****
      if ((-m1) + m2 - m3 .ne. 0) then
         gen_gaunt3_coeff_sub = 0.0d0
         return 
      endif

      !**** no triangle rule ****

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         gen_gaunt3_coeff_sub = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l2.eq.0 .and. l3.eq.0) then
         gen_gaunt3_coeff_sub = 0.0d0
         return 
      endif
      order = (l1 + l2 + l3) 

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL gauss_chebyshev_weights(x,w,order)
      gen_gaunt3_coeff_sub = 0.0d0
      do i = 1, order

        if ((m2.eq.0).and.(m3.eq.0)) then
          f = 0.0d0
        else if (m2.eq.0) then
         f = theta_lm(l1,m1,x(i))*
     >      (m3*theta_lm(l2,m2,x(i))* theta_lm(l3,m3,x(i)))
        else if (m3.eq.0) then
         f = theta_lm(l1,m1,x(i))*
     >      (m2*theta_lm(l2,m2,x(i))* theta_lm(l3,m3,x(i)))
        else
         f = theta_lm(l1,m1,x(i))*(m2+m3)*
     >      (theta_lm(l2,m2,x(i))* theta_lm(l3,m3,x(i)))
        end if

        gen_gaunt3_coeff_sub = gen_gaunt3_coeff_sub + w(i)*f

      end do

      gen_gaunt3_coeff_sub = gen_gaunt3_coeff_sub*2.0d0*pi

      return
      end 
 


*     *********************************************
*     *                                           *
*     *        paw_addto_potential_gntxc          *
*     *                                           *
*     *********************************************

      subroutine paw_init_gntxc(lmax)
      implicit none
      integer lmax

#include "paw_ma.fh"
#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      integer nkatm,nsize,nsize0
      integer nlmi,nilm,njlm,nb,nb2,lmax2
      integer ia,i,j,l,m,li,mi,lj,mj
      integer i_mtr,i_pot
      double precision gnt1

*     **** external functions ****
      double precision gen_gaunt_coeff
      external         gen_gaunt_coeff


      nkatm = ion_nkatm()
      lmax2 = (lmax+1)**2
      nsize = 0
      do ia=1,nkatm
        nb  = paw_basis_nbasis(ia)
        nb2 = nb*nb

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( ((l.le.(li+lj)) .and. (l.ge.ABS(li-lj))).or.
     >       (mod(l+li+lj,2).eq.0)                        ) then
         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt_coeff(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               nsize = nsize + 1
             end if
           end if
           end do !mi
           end do !mj
         end do !m
         end if
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l
      end do !ia


*     ***** allocate memory ****
      ok = .true.
      ok = ok .and. my_alloc(mt_int,nkatm,"ncount_gntxc",  ncount_gntxc)
      ok = ok .and. my_alloc(mt_int,nkatm,"indx_gntxc",indx_start_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"imtr_gntxc",      imtr_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"ipot_gntxc",      ipot_gntxc)
      ok = ok .and. my_alloc(mt_dbl,nsize,"ipot_gntxc",       gnt_gntxc)
      if (.not.ok)
     > call errquit("paw_init_gntxc: error allocating from heap",0,0)




*     **** xc potential  non-local matrix elements ****
      nsize = 0
      do ia=1,nkatm
         nb  = paw_basis_nbasis(ia)
         nb2 = nb*nb

         nsize0 = nsize
         int_mb(indx_start_gntxc(1)+ia-1) = nsize

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( ((l.le.(li+lj)) .and. (l.ge.ABS(li-lj))).or.
     >       (mod(l+li+lj,2).eq.0)                        ) then

         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt_coeff(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               i_mtr = (j-1)+(i-1)*nb
     >               + (nlmi+l+m)*nb2
               i_pot = (njlm+lj+mj)
     >               + (nilm+li+mi)*paw_proj_nbasis(ia)
               int_mb(imtr_gntxc(1)+nsize) = i_mtr+1
               int_mb(ipot_gntxc(1)+nsize) = i_pot+1
               dbl_mb(gnt_gntxc(1) +nsize)  = gnt1
               nsize = nsize + 1
             end if
 
           end if
           end do !mi
           end do !mj
         end do !m
 
         end if
 
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l

         int_mb(ncount_gntxc(1)+ia-1) = nsize-nsize0

      end do !ia

      return
      end

      subroutine paw_end_gntxc()
      implicit none

#include "paw_ma.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      ok = .true.
      ok = ok .and. my_dealloc(gnt_gntxc)
      ok = ok .and. my_dealloc(imtr_gntxc)
      ok = ok .and. my_dealloc(ipot_gntxc)
      ok = ok .and. my_dealloc(indx_start_gntxc)
      ok = ok .and. my_dealloc(ncount_gntxc)
      if (.not.ok)
     > call errquit("paw_end_gntxc: error freeing heap",0,0)

      return
      end


      subroutine paw_addto_potential_gntxc(ia,matr,v)
      implicit none
      integer ia
      double complex matr(*)
      double complex v(*)

#include "bafdecls.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >                   
*     **** local variables ****
      integer indx,ncount

      ncount = int_mb(ncount_gntxc(1)    +ia-1)
      indx   = int_mb(indx_start_gntxc(1)+ia-1)
      

      call paw_addto_potential_gntxc_sub(ncount,
     >                               dbl_mb(gnt_gntxc(1)+indx),
     >                               matr,int_mb(imtr_gntxc(1)+indx),
     >                               v,   int_mb(ipot_gntxc(1)+indx))

      return
      end

      subroutine paw_addto_potential_gntxc_sub(n,gnt,
     >                                        matr,indx_m,
     >                                        v,   indx_v)
      implicit none
      integer n
      double precision gnt(*)
      double complex matr(*)
      integer indx_m(*)
      double complex v(*)
      integer indx_v(*)

c     **** local varialbles ****
      integer i,ii,jj

      
      do i=1,n
         ii = indx_v(i)
         jj = indx_m(i)
         v(ii) = v(ii) + gnt(i)*matr(jj)
      end do
      return
      end







*     *********************************************
*     *                                           *
*     *        paw_addto_potential_gntxc2         *
*     *                                           *
*     *********************************************

      subroutine paw_init_gntxc2(lmax)
      implicit none
      integer lmax

#include "paw_ma.fh"
#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc2/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      integer nkatm,nsize,nsize0
      integer nlmi,nilm,njlm,nb,nb2,lmax2
      integer ia,i,j,l,m,li,mi,lj,mj
      integer i_mtr,i_pot
      double precision gnt1

*     **** external functions ****
      double precision gen_gaunt2_coeff_sub  !*** don't allocate gaunt tmp space
      external         gen_gaunt2_coeff_sub


      nkatm = ion_nkatm()
      lmax2 = (lmax+1)**2
      nsize = 0
      do ia=1,nkatm
        nb  = paw_basis_nbasis(ia)
        nb2 = nb*nb

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( (mod(l+li+lj,2).eq.1) ) then
         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt2_coeff_sub(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               nsize = nsize + 1
             end if
           end if
           end do !mi
           end do !mj
         end do !m
         end if
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l
      end do !ia


*     ***** allocate memory ****
      ok = .true.
      ok = ok .and. my_alloc(mt_int,nkatm,"ncount_gntxc",  ncount_gntxc)
      ok = ok .and. my_alloc(mt_int,nkatm,"indx_gntxc",indx_start_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"imtr_gntxc",      imtr_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"ipot_gntxc",      ipot_gntxc)
      ok = ok .and. my_alloc(mt_dbl,nsize,"ipot_gntxc",       gnt_gntxc)
      if (.not.ok)
     > call errquit("paw_init_gntxc2: error allocating from heap",0,0)




*     **** xc potential  non-local matrix elements ****
      nsize = 0
      do ia=1,nkatm
         nb  = paw_basis_nbasis(ia)
         nb2 = nb*nb

         nsize0 = nsize
         int_mb(indx_start_gntxc(1)+ia-1) = nsize

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( (mod(l+li+lj,2).eq.1) ) then

         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt2_coeff_sub(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               i_mtr = (j-1)+(i-1)*nb
     >               + (nlmi+l+m)*nb2
               i_pot = (njlm+lj+mj)
     >               + (nilm+li+mi)*paw_proj_nbasis(ia)
               int_mb(imtr_gntxc(1)+nsize) = i_mtr+1
               int_mb(ipot_gntxc(1)+nsize) = i_pot+1
               dbl_mb(gnt_gntxc(1) +nsize)  = gnt1
               nsize = nsize + 1
             end if
 
           end if
           end do !mi
           end do !mj
         end do !m
 
         end if
 
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l

         int_mb(ncount_gntxc(1)+ia-1) = nsize-nsize0

      end do !ia

      return
      end

      subroutine paw_end_gntxc2()
      implicit none

#include "paw_ma.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc2/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      ok = .true.
      ok = ok .and. my_dealloc(gnt_gntxc)
      ok = ok .and. my_dealloc(imtr_gntxc)
      ok = ok .and. my_dealloc(ipot_gntxc)
      ok = ok .and. my_dealloc(indx_start_gntxc)
      ok = ok .and. my_dealloc(ncount_gntxc)
      if (.not.ok)
     > call errquit("paw_end_gntxc2: error freeing heap",0,0)

      return
      end


      subroutine paw_addto_potential_gntxc2(ia,matr,v)
      implicit none
      integer ia
      double complex matr(*)
      double complex v(*)

#include "bafdecls.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc2/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >                   
*     **** local variables ****
      integer indx,ncount

      ncount = int_mb(ncount_gntxc(1)    +ia-1)
      indx   = int_mb(indx_start_gntxc(1)+ia-1)
      

      call paw_addto_potential_gntxc2_sub(ncount,
     >                               dbl_mb(gnt_gntxc(1)+indx),
     >                               matr,int_mb(imtr_gntxc(1)+indx),
     >                               v,   int_mb(ipot_gntxc(1)+indx))

      return
      end

      subroutine paw_addto_potential_gntxc2_sub(n,gnt,
     >                                        matr,indx_m,
     >                                        v,   indx_v)
      implicit none
      integer n
      double precision gnt(*)
      double complex matr(*)
      integer indx_m(*)
      double complex v(*)
      integer indx_v(*)

c     **** local varialbles ****
      integer i,ii,jj

      
      do i=1,n
         ii = indx_v(i)
         jj = indx_m(i)
         v(ii) = v(ii) + gnt(i)*matr(jj)
      end do
      return
      end










*     *********************************************
*     *                                           *
*     *        paw_addto_potential_gntxc3         *
*     *                                           *
*     *********************************************

      subroutine paw_init_gntxc3(lmax)
      implicit none
      integer lmax

#include "paw_ma.fh"
#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc3/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      integer nkatm,nsize,nsize0
      integer nlmi,nilm,njlm,nb,nb2,lmax2
      integer ia,i,j,l,m,li,mi,lj,mj
      integer i_mtr,i_pot
      double precision gnt1

*     **** external functions ****
      double precision gen_gaunt3_coeff_sub  !*** don't allocate gaunt tmp space
      external         gen_gaunt3_coeff_sub


      nkatm = ion_nkatm()
      lmax2 = (lmax+1)**2
      nsize = 0
      do ia=1,nkatm
        nb  = paw_basis_nbasis(ia)
        nb2 = nb*nb

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( (mod(l+li+lj,2).eq.0) ) then
         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt3_coeff_sub(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               nsize = nsize + 1
             end if
           end if
           end do !mi
           end do !mj
         end do !m
         end if
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l
      end do !ia


*     ***** allocate memory ****
      ok = .true.
      ok = ok .and. my_alloc(mt_int,nkatm,"ncount_gntxc",  ncount_gntxc)
      ok = ok .and. my_alloc(mt_int,nkatm,"indx_gntxc",indx_start_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"imtr_gntxc",      imtr_gntxc)
      ok = ok .and. my_alloc(mt_int,nsize,"ipot_gntxc",      ipot_gntxc)
      ok = ok .and. my_alloc(mt_dbl,nsize,"ipot_gntxc",       gnt_gntxc)
      if (.not.ok)
     > call errquit("paw_init_gntxc3: error allocating from heap",0,0)




*     **** xc potential  non-local matrix elements ****
      nsize = 0
      do ia=1,nkatm
         nb  = paw_basis_nbasis(ia)
         nb2 = nb*nb

         nsize0 = nsize
         int_mb(indx_start_gntxc(1)+ia-1) = nsize

         nlmi = 0
         do l=0,lmax
         nilm = 0
         do i=1,nb
         li = paw_basis_orb_l(i,ia)
         njlm = 0
         do j=1,nb
         lj = paw_basis_orb_l(j,ia)
         if( (mod(l+li+lj,2).eq.0) ) then

         do m=-l,l
           do mi=-li,li
           do mj=-lj,lj
            if(m.eq.(mi-mj)) then
             gnt1 = gen_gaunt3_coeff_sub(l,m,li,mi,lj,mj)
             if (abs(gnt1).gt.1.0d-9) then
               i_mtr = (j-1)+(i-1)*nb
     >               + (nlmi+l+m)*nb2
               i_pot = (njlm+lj+mj)
     >               + (nilm+li+mi)*paw_proj_nbasis(ia)
               int_mb(imtr_gntxc(1)+nsize) = i_mtr+1
               int_mb(ipot_gntxc(1)+nsize) = i_pot+1
               dbl_mb(gnt_gntxc(1) +nsize)  = gnt1
               nsize = nsize + 1
             end if
 
           end if
           end do !mi
           end do !mj
         end do !m
 
         end if
 
         njlm = njlm + 2*lj+1
         end do !j
         nilm = nilm + 2*li+1
         end do !i
         nlmi=nlmi+2*l+1
         end do !l

         int_mb(ncount_gntxc(1)+ia-1) = nsize-nsize0

      end do !ia

      return
      end

      subroutine paw_end_gntxc3()
      implicit none

#include "paw_ma.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc3/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >

      logical ok
      ok = .true.
      ok = ok .and. my_dealloc(gnt_gntxc)
      ok = ok .and. my_dealloc(imtr_gntxc)
      ok = ok .and. my_dealloc(ipot_gntxc)
      ok = ok .and. my_dealloc(indx_start_gntxc)
      ok = ok .and. my_dealloc(ncount_gntxc)
      if (.not.ok)
     > call errquit("paw_end_gntxc3: error freeing heap",0,0)

      return
      end


      subroutine paw_addto_potential_gntxc3(ia,matr,v)
      implicit none
      integer ia
      double complex matr(*)
      double complex v(*)

#include "bafdecls.fh"
#include "errquit.fh"

      integer gnt_gntxc(4)
      integer imtr_gntxc(4)
      integer ipot_gntxc(4)
      integer indx_start_gntxc(4)
      integer ncount_gntxc(4)
      common /paw_gntxc3/ gnt_gntxc,imtr_gntxc,ipot_gntxc,
     >                    indx_start_gntxc,ncount_gntxc
     >                   
*     **** local variables ****
      integer indx,ncount

      ncount = int_mb(ncount_gntxc(1)    +ia-1)
      indx   = int_mb(indx_start_gntxc(1)+ia-1)
      

      call paw_addto_potential_gntxc3_sub(ncount,
     >                               dbl_mb(gnt_gntxc(1)+indx),
     >                               matr,int_mb(imtr_gntxc(1)+indx),
     >                               v,   int_mb(ipot_gntxc(1)+indx))

      return
      end

      subroutine paw_addto_potential_gntxc3_sub(n,gnt,
     >                                        matr,indx_m,
     >                                        v,   indx_v)
      implicit none
      integer n
      double precision gnt(*)
      double complex matr(*)
      integer indx_m(*)
      double complex v(*)
      integer indx_v(*)

c     **** local varialbles ****
      integer i,ii,jj

      
      do i=1,n
         ii = indx_v(i)
         jj = indx_m(i)
         v(ii) = v(ii) + dcmplx(0.0d0,gnt(i))*matr(jj)
      end do
      return
      end
