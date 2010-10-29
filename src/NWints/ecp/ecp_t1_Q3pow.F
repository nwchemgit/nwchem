C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q3pow (n,m,x,xt,gamma,term,test,Q,tol)
*                                                                      *
*   Calculate Q^{m+p}_m for p = 3 using the series expansion.          *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q functions to be evaluated for given order    *
*   m (inp) - order of Q function                                      *
*   x (inp) - array of values beta/2sqrt(alpha)                        *
*   xt (scr) - array to store x^2                                      *
*   gamma (scr) - array of values of 1/sqrt(alpha)                     *
*   term (scr) - array of terms in series                              *
*   test (scr) - array of ratios of terms to sums                      *
*   Q (out) - Q^m+3_m functions                                        *
*   tol (inp) - maximum relative error in Q functions                  *
*                                                                      *
*   Written by W. A. de Jong                                           *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
#include "errquit.fh"
      integer i,idamax,m,n
      double precision x(n),xt(n),gamma(n),term(n),test(n),Q(n),
     &    tol,fac0,fac1,fac2,fac3
*
      if (m .lt. 0) call errquit ('Negative m in ecp_t1_Q3pow',99,
     &       BASIS_ERR)
      fac3 = half
      fac1 = zero
      fac2 = one
      do i = 1,m
        fac1 = fac1+two
        fac2 = fac2+two  
        fac3 = fac3*fac1/fac2
      end do
      fac3 = fac3*(m+1)
      if (m .eq. 0) then
        do i = 1,n
          Q(i) = fac3*(gamma(i)**4)
        end do
      else
        do i = 1,n
          Q(i) = fac3*(gamma(i)**4)*(x(i)*gamma(i))**m
        end do
      end if
      fac1 = fac1+four 
      fac2 = fac2+two              
      fac0 = fac1/fac2
      do i = 1,n
        xt(i) = x(i)*x(i)
        test(i) = fac0*xt(i)
        term(i) = Q(i)*test(i)
      end do
      fac3 = 1
    1 i = idamax(n,test,1)
      if (test(i) .ge. tol) then
        fac1 = fac1+two
        fac2 = fac2+two
        fac3 = fac3+one
        fac0 = fac1/(fac2*fac3)
        do i = 1,n
          if (test(i) .ge. tol) then
            Q(i) = Q(i)+term(i)
            term(i) = term(i)*xt(i)*fac0
            test(i) = term(i)/Q(i)
          end if
        end do
        go to 1
      end if
*
      do i = 1,n
        Q(i) = Q(i)*exp(-xt(i))
      end do
*
      return
      end
