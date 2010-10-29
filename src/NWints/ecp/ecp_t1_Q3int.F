C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q3int (n,m,x,a_sqrt,Q3,temp,ind,tol)
*                                                                      *
*   Driver for calculation of type 1 Q^{m+p}_{m} integrals for p = 3   *
*                                                                      *
* Q^{m+p}_{m} = \int_0^\infty dr r^{m+p} e^{-\alpha r^2} i_m(\beta r)  *
*                                                                      *
*   where $i_m(\beta r)$ is the modified spherical Bessel function.    *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n (inp) - number of Q integrals to be evaluated for given order    *
*   m (inp) - order of Q integral                                      *
*   x (inp) - array of values beta/2sqrt(alpha)                        *
*   a_sqrt (scr) - array of values of one/sqrt(alpha)                  *
*   Q3 (out) - Q^{m}_{m} integrals                                     *
*   temp (scr) - scratch array                                         *
*   ind (scr) - index array                                            *
*   tol (inp) - maximum relative error in bessel functions             *
*                                                                      *
*   Written by W. A. de Jong                                           *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,m,n,na,np,n1,ind(n)
      double precision x(n),a_sqrt(n),Q3(n),temp(n,6),tol,
     &    m_plus_half,m_plus_2,log_accy
*
*   Initialize ind array
*
      do i = 1,n
        ind(i) = -1
      enddo
*
*   Gather arguments
*
      m_plus_2 = m+2
      m_plus_half = m_plus_2-half
      log_accy = log(tol)-log(two)-m_plus_half*log(m_plus_2)+m_plus_2
      na = 0
      np = 0
      n1 = n+1
      do i = 1,n
        if (x(i) .eq. zero) then
          if (m .eq. 0) then
            Q3(i) = a_sqrt(i)**4/two
          else
            Q3(i) = zero
          end if
        else if (log_accy+x(i)**2+two*m_plus_half*log(x(i)) .gt. zero) 
     &        then
          na = na+1
          temp(na,1) = x(i)
          temp(na,3) = a_sqrt(i)
          ind(na) = i
        else
          np = np+1
          temp(n1-np,1) = x(i)
          temp(n1-np,3) = a_sqrt(i)
          ind(n1-np) = i
        end if
      end do
      i = n1-np
*
      if (na.gt.0.or.np.gt.0) then
        if (na .gt. 0) call ecp_t1_Q3asy (na,m,temp(1,1),temp(1,2),
     &      temp(1,3),temp(1,6),tol)
        if (np .gt. 0) call ecp_t1_Q3pow (np,m,temp(i,1),temp(i,2),
     &      temp(i,3),temp(i,4),temp(i,5),temp(i,6),tol)
        do i = 1,n
          if (ind(i).gt.0)
     &        Q3(ind(i)) = temp(i,6)
        end do
      endif
*
      return
      end
