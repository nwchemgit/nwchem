c
c $Id: paw_gaunt.F,v 1.4 2003-10-21 02:05:16 marat Exp $
c
c $Log: not supported by cvs2svn $
c Revision 1.3  2003/02/11 01:41:23  edo
c eliminated f90-isms
c
c Revision 1.2  2003/02/10 02:36:26  marat
c ...
c
c Revision 1.1  2003/02/06 06:04:37  marat
c Initial Revision
c
      subroutine paw_gaunt_init()
      !
      implicit none
#include "paw_gaunt_data.fh"

      paw_gaunt_lmax = 4
      
      return
      end
      function gen_gaunt_coeff (l1, m1, l2, m2, l3, m3)

      implicit none
      integer l1, l2, l3
      integer m1, m2, m3
      integer i
      integer order
      double precision x1, x2
      double precision  gen_gaunt_coeff

#include"paw_params.fh"
#include"paw_spher_func.fh"

c     work arrays for integration
      double precision x(100),w(100)

      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) call errquit(
     1   'Invalid parameter in gen_gaunt_coeff, negative l',0,1)
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2)) call 
     1   errquit('Invalid parameter in gen_gaunt_coeff, m > l',0,1)
!    Do integration over angle phi
      if ((-m1) + m2 - m3 .ne. 0) then
         gen_gaunt_coeff = 0.0
         return 
      endif
!    Check the triangle rule
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         gen_gaunt_coeff = 0.0
         return 
      endif
!    Check if the integrand is odd function==>integral is zero
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         gen_gaunt_coeff = 0.0
         return 
      endif
!    hANDLE THE EXEPTIONAL CASE
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         gen_gaunt_coeff = 1.0/sqrt(fourpi)
         return 
      endif
      x1 = -1.0
      x2 = 1.0
      order = l1 + l2 + l3
!    Generate weights and coordinates for Gauss-Legendre integration
      CALL gauss_weights(x1,x2,x,w,order)
      gen_gaunt_coeff = 0.0
      do i = 1, order
         gen_gaunt_coeff = gen_gaunt_coeff + w(i)*theta_lm(l1,m1,x(i))*
     1      theta_lm(l2,m2,x(i))*theta_lm(l3,m3,x(i))
      end do
      gen_gaunt_coeff = gen_gaunt_coeff*2.0*pi
      end 
 
