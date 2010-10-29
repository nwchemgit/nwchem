C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_down_m (m_min,m_max,n,a_inv,Q_hi_a,Q_mid_a,Q_lo_a)
*                                                                      *
*   Perform downward recursion in m using the relation for bessel      *
*   functions i_{m-1}(ar) = i_{m+1}(ar) + (2m+1)/(ar) i_{m}(ar)        *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   m_min (inp) - minimum value of m                                   *
*   m_max (inp) - maximum value of m                                   *
*   n (inp) - number of functions                                      *
*   a_inv (inp) - inverse of factors multiplying r                     *
*   Q_hi_a (inp) - integral containing i_{m+1}                         *
*   Q_mid_a (inp) - integral containing i_{m}/r                        *
*   Q_lo_a (out) - integral containing i_{m-1}                         *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,m,m_min,m_max,n
      double precision a_inv(n),Q_hi_a(n,m_min:m_max),
     &    Q_mid_a(n,m_min:m_max),Q_lo_a(n,m_min:m_max),fac
*
      do m = m_max,m_min,-1
        fac = 2*m+1
        do i = 1,n
          Q_lo_a(i,m) = Q_hi_a(i,m) + fac*a_inv(i)*Q_mid_a(i,m)
        end do
      end do
*
      return
      end
