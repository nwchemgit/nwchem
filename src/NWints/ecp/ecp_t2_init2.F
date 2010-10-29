C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t2_init2 (n,ldQ,m_min,m_max,h,
     &    tol,alpha,beta,gamma,prefactor,temp,x,ind,Q)
*                                                                      *
*   Set up values of Q^2_mm, Q^3_m-1m and Q^3_mm-1 to initiate the     *
*   recursion                                                          *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of functions to evaluate                          *
*   ldQ (inp) - leading dimension of array Q                           *
*   m_min (inp) - minimum value of m                                   *
*   m_max (inp) - maximum value of m                                   *
*   h (inp) switch between alpha and beta for l_b < l_a                *
*   tol (inp) - maximum relative error in bessel functions             *
*   alpha (inp) -  a/2sqrt(c)                                          *
*   beta (inp) - b/2sqrt(c)                                            *
*   gamma (inp) - 1/sqrt(c)                                            *
*   prefactor (inp) - exponential prefactor (see calling routine)      *
*   temp - work array                                                  *
*   x (scr) - array of values ab/2c                                    *
*   ind (scr) - index array for bessel driver                          *
*   Q (out) - uncontracted Q integrals                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer h,i,ldQ,m,m_max,m_min,n,ind(n)
      double precision alpha(n),beta(n),gamma(n),temp(n*5),x(n),
     &    prefactor(n),Q(ldQ,m_min:m_max,2),tol,fac
*
*   Set up argument values
*
      if (n .eq. 0) return
      do i = 1,n
        x(i) = two*alpha(i)*beta(i)
      end do
*
*   Evaluate bessel functions for Q^2_{mm} with m = m_max, m_max-1
*
      do m = m_max,max(m_max-1,m_min),-1
        call ecp_bessel (n,m,x,Q(1,m,1),temp,ind,tol)
*
*   Multiply by (sqrt{pi}/4) c^{-3/2} and prefactor
*
        do i = 1,n
          Q(i,m,1) = qpirt*prefactor(i)*Q(i,m,1)*gamma(i)**3
        end do
      end do
*
*   Downward recursion for lower Q^2.
*
      do m = m_max-2, m_min, -1
        fac = m+m+3
        do i = 1,n
          Q(i,m,1) = Q(i,m+2,1)+fac*Q(i,m+1,1)/x(i)
        end do
      end do
      if (m_min .eq. m_max) return
*
*   Odd k recursion to obtain Q^3_{m+1m} or Q^3_{mm+1}
*
      call ecp_up_k (m_min,m_max,1,1,h,n,ldQ,alpha,beta,gamma,
     &    Q(1,m_min,1),Q(1,m_min,1),Q(1,m_min,2))
*
      return
      end
