C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q0pow (n,m,x,xt,a_sqrt,term,test,Q,tol)
*                                                                      *
*   Calculate Q^{m+p}_m for p = 0 using the series expansion.          *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q functions to be evaluated for given order    *
*   m (inp) - order of Q function                                      *
*   x (inp) - array of values beta/2sqrt(alpha)                        *
*   xt (scr) - array to store x^2                                      *
*   a_sqrt (scr) - array of values of 1/sqrt(alpha)                    *
*   term (scr) - array of terms in series                              *
*   test (scr) - array of ratios of terms to sums                      *
*   Q (out) - Q^m_m functions                                          *
*   tol (inp) - maximum relative error in Q functions                  *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
#include "errquit.fh"
      integer i,j,idamax,n,m
      double precision x(n),xt(n),a_sqrt(n),term(n),test(n),Q(n),
     &    tol,fac0,fac1,fac2,fac3
*
      if (m .lt. 0) call errquit ('Negative m in ecp_t1_Q0pow',99,
     &       BASIS_ERR)
      do i = 1,n
        xt(i) = x(i)*x(i)
        Q(i) = zero
        test(i) = one
      end do
      if (m .eq. 0) then
        do i = 1,n
          term(i) = qpirt*a_sqrt(i)
        end do
      else
        do i = 1,n
          term(i) = qpirt*a_sqrt(i)*(x(i)*a_sqrt(i))**m
        end do
      end if
      fac1 = m+m+1
      fac1 = fac1/two
      fac2 = zero
    1 j = idamax(n,test,1)
      if (test(j) .ge. tol) then
        fac0 = one/fac1
        fac1 = fac1+one
        fac2 = fac2+one
        fac3 = one/fac2
        do i = 1,n
          if (test(i) .ge. tol) then
            Q(i) = Q(i)+term(i)*fac0
            term(i) = term(i)*xt(i)*fac3
            test(i) = term(i)/(Q(i)*fac1)
          end if
        end do
        go to 1
      end if
C      j = fac2+0.5d0
C      write (6,*) ' Terms in Q^m_m power series: ',j
*
      do i = 1,n
        Q(i) = Q(i)*exp(-xt(i))
      end do
*
      return
      end
