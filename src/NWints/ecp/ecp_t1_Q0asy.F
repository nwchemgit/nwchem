C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q0asy (n,m,x,xt,a_sqrt,term,test,Q,tol)
*                                                                      *
*   Calculate Q^{m+p}_m function for p = 0 by use of the asymptotic    * 
*   expansion                                                          *
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
*   tol (inp) - maximum relative error in Q0 functions                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,idamax,m,n
      double precision x(n),xt(n),a_sqrt(n),term(n),test(n),Q(n),tol,
     &    fac1
*
      do i = 1,n
        xt(i) = x(i)*x(i)
        term(i) = qpirt*a_sqrt(i)*(x(i)*a_sqrt(i))**m/xt(i)
        Q(i) = zero
        test(i) = one
      end do
      fac1 = 1-m-m
      fac1 = fac1/two
*
    1 i = idamax(n,test,1)
      if (test(i) .ge. tol) then
CDIR$ NOVECTOR
        do i = 1,n
          if (test(i) .ge. tol) then
            Q(i) = Q(i)+term(i)
            term(i) = (term(i)*fac1)/xt(i)
            test(i) = sign(abs(term(i))/Q(i),xt(i)-fac1)
          end if
        end do
        fac1 = fac1+one
        go to 1
      end if
C      i = fac1+m
C      write (6,*) ' Terms in Q^m_m asymptotic series: ',i
*
      return
      end
