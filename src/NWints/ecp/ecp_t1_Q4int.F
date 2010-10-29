C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q4int (n,m,x,a_sqrt,Q4)
*                                                                      *
*   Routine to evaluate type 1 Q^{m+p}_{m} integrals for p = 4         *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q integrals to be evaluated for given order    *
*   m (inp) - order of Q integral                                      *
*   x (inp) - array of values beta/2sqrt(alpha)                        *
*   a_sqrt (scr) - array of values of one/sqrt(alpha)                  *
*   Q4 (out) - Q^{m+4}_{m} integrals                                   *
*                                                                      *
*   Written by W. A. de Jong                                          *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,n,m
      double precision x(n),a_sqrt(n),Q4(n)
*
      if (m .eq. 0) then
        do i = 1,n
          Q4(i) = qpirt*(a_sqrt(i)**5)*((x(i)**2)+1.5d0)
        end do
      else
        do i = 1,n
          Q4(i) = qpirt*(a_sqrt(i)**(m+5))
          Q4(i) = Q4(i)*((x(i)**(m+2))+(x(i)**m)*(m+1.5d0))
        end do
      end if
*
      return
      end
