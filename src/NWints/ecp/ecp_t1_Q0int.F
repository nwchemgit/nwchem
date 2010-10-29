C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t1_Q0int (n,m,x,a_sqrt,Q0,temp,ind,tol)
*                                                                      *
*   Driver for calculation of type 1 Q^{m+p}_{m} integrals for p = 0   *
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
*   a_sqrt (scr) - array of values of 1/sqrt(alpha) for each integral  *
*   Q0 (out) - Q^{m}_{m} integrals                                     *
*   temp (scr) - scratch array                                         *
*   ind (scr) - index array                                            *
*   tol (inp) - maximum relative error in bessel functions             *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer i,m,n,na,np,n1,ind(n)
      double precision x(n),a_sqrt(n),q0(n),temp(n,6),tol,
     &    m_plus_half,real_m,log_accy
*
*   Initialize ind array
*
      do i = 1,n
        ind(i) = -1
      enddo
*
*   Gather arguments
*
      real_m = m
      m_plus_half = real_m+half
      log_accy = log(tol)-log(two)-real_m*log(m_plus_half)+m_plus_half
      na = 0
      np = 0
      n1 = n+1
      do i = 1,n
        if (x(i) .eq. zero) then
          if (m .eq. 0) then
            Q0(i) = qpirt*two*a_sqrt(i)
          else
            Q0(i) = zero
          end if
        else if (log_accy+x(i)**2+two*real_m*log(x(i)) .gt. zero) then
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
*debug      if (na.gt.0.or.np.gt.0) then
*debug        write (6,*) 'Q0:na,np,sum,n,diff',na,np,(na+np),n,(n-na-np)
*debug      endif
*
      if (na.gt.0.or.np.gt.0) then
*................... i.e., ind(*) is actually set by conditions above !
        if (na .gt. 0) call ecp_t1_Q0asy (na,m,temp(1,1),temp(1,2),
     &      temp(1,3),temp(1,4),temp(1,5),temp(1,6),tol)
        if (np .gt. 0) call ecp_t1_Q0pow (np,m,temp(i,1),temp(i,2),
     &      temp(i,3),temp(i,4),temp(i,5),temp(i,6),tol)
        do i = 1,n
          if (ind(i).gt.0)
     &        Q0(ind(i)) = temp(i,6)
        end do
      endif
*
      return
      end
