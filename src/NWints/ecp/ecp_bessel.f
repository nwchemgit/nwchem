C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_bessel (n,m,x,bessel,temp,ind,tol)
*                                                                      *
*   Calculate modified spherical bessel function exp(-x) i_m(x)        *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of bessel functions to be evaluated for given     *
*             order                                                    *
*   m (inp) - order of bessel function                                 *
*   x (inp) - array of arguments for bessel function                   *
*   bessel (out) - bessel functions                                    *
*   temp (scr) - scratch array                                         *
*   ind (scr) - index array                                            *
*   tol (inp) - maximum relative error in bessel functions             *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer n,m,i,na,np,ind(n)
      double precision x(n),bessel(n),temp(n,5),tol,two,cut0,cutoff
      parameter (two = 2.0d00, cut0 = 0.106d0)
*
*   Determine cutoff for division between power and asymptotic series
*
      if (m .eq. 0) then
        cutoff = cut0
      else
        cutoff = m*(m+1)/2
      end if
      cutoff = cutoff/two
*
*   Gather arguments for power and asymptotic series
*
      na = 0
      np = n+1
      do i = 1,n
        if (x(i) .gt. cutoff) then
          na = na+1
          temp(na,1) = x(i)
          ind(na) = i
        else
          np = np-1
          temp(np,1) = x(i)
          ind(np) = i
        end if
      end do
      i = np
      np = n-na
*
*   Evaluate functions and scatter into output array
*
      if (np .gt. 0) call ecp_bess_pow (np,m,temp(i,1),temp(i,2),
     &    temp(i,3),temp(i,4),temp(i,5),tol)
      if (na .gt. 0) call ecp_bess_asy (na,m,temp(1,1),temp(1,2),
     &    temp(1,3),temp(1,4))
      do i = 1,n
        bessel(ind(i)) = temp(i,2)
      end do
*
      return
      end
