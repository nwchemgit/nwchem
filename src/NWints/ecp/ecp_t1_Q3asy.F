C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q3asy (n,m,x,xt,gamma,Q,tol)
*                                                                      *
*   Calculate Q^{m+p}_m function for p = 3 by use of the asymptotic    * 
*   expansion                                                          *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q functions to be evaluated for given order    *
*   m (inp) - order of Q function                                      *
*   x (inp) - array of values beta/(2*sqrt(alpha))                     *
*   xt (scr) - array to store x^2                                      *
*   gamma (scr) - array of values of 1/sqrt(alpha)                     *
*   Q (out) - Q^m_m functions                                          *
*   tol (inp) - maximum relative error in Q3 functions                 *
*                                                                      *
*   Written by W. A. de Jong                                           *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,k,m,n
      double precision x(n),xt(n),gamma(n),Q(n),tol,fac0,fac1,fac2,fac3
*
      do i = 1,n
        xt(i) = one/(x(i)*x(i))
        Q(i) = one
      end do
*
      fac1 = half*(2*m-1)
      fac2 = -m
      fac3 = m+1
      do k = 0,m-1
        fac0 = fac1*fac2/fac3
        do i = 1,n
          Q(i) = one+Q(i)*fac0*xt(i)
        end do
        fac1 = fac1 - one
        fac2 = fac2 - one
        fac3 = fac3 - one
      end do
      do i = 1,n
        Q(i) = Q(i)*qpirt*(x(i)*gamma(i))**(m+1)*gamma(i)**3
      end do
*
      return
      end
