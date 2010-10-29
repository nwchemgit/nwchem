C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q1asy (n,m,x,xt,gamma,Q,tol)
*                                                                      *
*   Calculate Q^{m+p}_m function for p = 1 by use of the asymptotic    * 
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
*   tol (inp) - maximum relative error in Q1 functions                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,k,m,n
      double precision x(n),xt(n),gamma(n),Q(n),tol,fac0,fac1,fac2
*
      do i = 1,n
        xt(i) = one/(x(i)*x(i))
        Q(i) = one
      end do
*
      do k = 0,m-1
        fac1 = k+1
        fac2 = 2*(m-k)
        fac0 = (fac2-one)*fac1/fac2
        do i = 1,n
          Q(i) = one-Q(i)*fac0*xt(i)
        end do
      end do
      do i = 1,n
        Q(i) = Q(i)*qpirt*(x(i)*gamma(i))**(m+1)*xt(i)*gamma(i)
      end do
*
      return
      end
