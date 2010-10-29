C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q2int (n,m,x,a_sqrt,Q2)
*                                                                      *
*   Routine to evaluate type 1 Q^{m+p}_{m} integrals for p = 2         *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q integrals to be evaluated for given order    *
*   m (inp) - order of Q integral                                      *
*   x (inp) - array of values beta/2sqrt(alpha)                        *
*   a_sqrt (scr) - array of values of one/sqrt(alpha)                  *
*   Q2 (out) - Q^{m+2}_{m} integrals                                   *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,n,m
      double precision x(n),a_sqrt(n),Q2(n)
*
      if (m .eq. 0) then
        do i = 1,n
          Q2(i) = qpirt*(a_sqrt(i)**3)
        end do
      else
        do i = 1,n
          Q2(i) = qpirt*(x(i)**m)*(a_sqrt(i)**(m+3))
        end do
      end if
*
      return
      end
