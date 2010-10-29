C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_bess_pow (n,m,x,bess,term,xt,test,tol)
*                                                                      *
*   Calculate modified spherical bessel function i_m(x) using the      *
*   series expansion.                                                  *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of bessel functions to be evaluated for given     *
*             order                                                    *
*   m (inp) - order of bessel function                                 *
*   x (inp) - array of arguments for bessel function                   *
*   bess (out) - bessel functions                                      *
*   term (scr) - array of terms in series                              *
*   xt (scr) - work space for functions of x                           *
*   test (scr) - array of ratios of terms to sums                      *
*   prefactor (scr) - x^m/(2m+1)!!                                     *
*   tol (inp) - maximum relative error in bessel functions             *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer n,m,i,idamax
      double precision x(n),bess(n),term(n),xt(n),test(n),
     &    tol,fac0,fac1,fac2,one,two
      parameter (one = 1.0d00, two = 2.0d00)
*
      if (n .le. 0) return
*
      do i = 1,n
        xt(i) = x(i)*x(i)
        bess(i) = 0.0d00
        term(i) = exp(-x(i))
        test(i) = one
      end do
      fac2 = 0.0d00
      fac1 = m+m+1
*
    1 fac1 = fac1+two
      fac2 = fac2+two
      fac0 = one/(fac1*fac2)
      do i = 1,n
        if (test(i) .ge. tol) then
          bess(i) = bess(i)+term(i)
          term(i) = term(i)*xt(i)*fac0
          test(i) = term(i)/bess(i)
        end if
      end do
      i = idamax(n,test,1)
      if (i .gt. 0) then
        if (test(i) .ge. tol) go to 1
      end if
      if (m .eq. 0) return
*
      fac0 = one
      do i = 3,m+m+1,2
        fac1 = i
        fac0 = fac0/fac1
      end do
      do i = 1,n
        bess(i) = bess(i)*fac0*x(i)**m
      end do
      return
      end
