C$Id: ecp_t2_p0pow.F,v 1.1 1996-10-16 22:52:31 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine ecp_t2_p0pow (n,m,k,alpha,beta,gamma,term,test,temp,
     &    Q,tol)
*                                                                      *
*   Calculate Q^|k|_{m+km} functions for k = 1,0,-1, and alpha < beta  *
*   by single power series expansion of the bessel function for alpha  *
*   in the radial integral.                                            *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q functions to be evaluated for given          *
*             order                                                    *
*   m (inp) - order of Q function                                      *
*   k (inp) - rank of Q function                                       *
*   alpha (inp) - array of values a/2sqrt(c)                           *
*   beta (inp) - array of values b/2sqrt(c)                            *
*   c - array of c values                                              *
*   term (scr) - array of terms alpha^j/(2j)!!(2m+2j+1)!!              *
*   test (scr) - array of ratios of terms to sums                      *
*   temp (scr) - work array, used for Q_m^{m+2j} recursion             *
*   Q (out) - Q functions                                              *
*   tol (inp) - maximum relative error in Q functions                  *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,idamax,j1,j2,j3,k,k_abs,m,n
      double precision alpha(n),beta(n),gamma(n),term(n),test(n),
     &    temp(n,9),Q(n),tol,fac0,fac1,fac2,fac3,fac4,fac5
*
*   Set up factors for recursion
*
      if (n .eq. 0) return
      fac1 = m+1+min(k,0)
      fac1 = fac1+half
      fac2 = 1-max(k,0)
      fac3 = zero
      fac4 = m
      fac4 = fac4+half
      fac5 = fac4+one
C      write (6,'(5F10.4)') fac1,fac2,fac3,fac4,fac5
*
*   First term for k < 1 is Q_m^m, which is not connected with the rest.
*   Series starts with Q_m^m+2
*
      do i = 1,n
        temp(i,1) = one
        Q(i) = zero
      end do
      if (k .le. 0) call ecp_t1_Q0int (n,m,beta,temp,Q,
     &      temp(1,3),temp(1,9),tol)
      call ecp_t1_Q2int (n,m,beta,temp,temp(1,2))
      fac0 = one/fac1
      do i = 1,n
        temp(i,4) = alpha(i)**2
        temp(i,5) = beta(i)**2
        if (k .eq. 1) then
          term(i) = one
        else
          term(i) = fac0*temp(i,4)
        end if
        test(i) = one
      end do
      j1 = 0
C      write (6,'(5F10.4)') fac1,fac2,fac3,fac4,fac5
    1 i = idamax(n,test,1)
      if (test(i) .ge .tol) then
        j1 = mod(j1,3)+1
        j2 = mod(j1,3)+1
        j3 = mod(j2,3)+1
        fac1 = fac1+one
        fac2 = fac2+one
        fac0 = one/(fac2*fac1)
        do i = 1,n
          if (test(i) .gt. tol) then
            Q(i) = Q(i)+temp(i,j2)*term(i)
            temp(i,j3) = (fac5+temp(i,5))*temp(i,j2)
     &                 - fac3*fac4*temp(i,j1)
            term(i) = term(i)*temp(i,4)*fac0
            test(i) = term(i)*temp(i,j3)/Q(i)
          end if
        end do
        fac3 = fac3+one
        fac4 = fac4+one
        fac5 = fac5+two
        go to 1
      end if
      i = fac3+0.1d0
C      write (6,*) ' Terms in single power series: ',i
C      write (6,*) 'beta,Q^0_mm values'
C      write (6,'(1p2e20.10)') (beta(i),Q(i),i=1,n)
*
*   Prefactor is (2*alpha)^{m-|k|}/(2(m-|k|)+1)!! /sqrt(c)^(m+|k|)
*
      k_abs = abs(k)
      fac1 = half
      fac5 = one
      do i = 1,m+k
        fac1 = fac1+one
        fac5 = fac5/fac1
      end do
C      write (6,*) fac1,fac5
      do i = 1,n
        Q(i) = Q(i)*alpha(i)**(m+k)*fac5*(gamma(i)**(k_abs+1))
      end do
*
*   Perform remaining exponent shift.
*
      do i = 1,n
        Q(i) = Q(i)*exp(-alpha(i)*(two*beta(i)+alpha(i)))
      end do
*
      return
      end
