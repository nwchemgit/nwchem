C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_bess_asy (n,m,x,bess,temp,xt)
*                                                                      *
*   Calculate modified spherical bessel function exp(-x) i_m(x) using  *
*   the asymptotic expansion.                                          *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of bessel functions to be evaluated for given     *
*             order                                                    *
*   m (inp) - order of bessel function                                 *
*   x (inp) - array of arguments for bessel function                   *
*   bess (out) - bessel functions                                      *
*   temp (scr) - work space for positive series                        *
*   xt (scr) - work space for functions of x                           *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
*
      integer i,n,m,j,mm
      double precision x(n),bess(n),temp(n),xt(n),top,bot,fact
*
      if (n .le. 0) return
*
      do i = 1,n
        bess(i) = one
        temp(i) = one
        xt(i) = one/x(i)
      end do
      fact = one
      mm = m+m+1
      do j = 1,m
        top = (mm-j)*j
        bot = 2*(m+1-j)
        fact = top/bot
        do i = 1,n
          bess(i) = one-fact*bess(i)*xt(i)
          temp(i) = one+fact*temp(i)*xt(i)
        end do
      end do
      fact = (-1)**(m+1)
      do i = 1,n
        xt(i) = x(i)+x(i)
        if (xt(i).gt.ln_thr_ecp) then  ! check for possible underflow problems
          bess(i) = (bess(i)+fact*temp(i)*exp(-xt(i)))/xt(i)
        else
          bess(i) = bess(i)/xt(i)
        endif
      end do
*
      return
      end
