C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_t2_ghq (k,m,n,nt,order,alpha,beta,gamma,temp,ind,
     &    Q,tol)
*                                                                      *
*   Evaluate modified Q^k_{mn} functions by Gauss-Hermite quadrature   *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   k (inp) - power of r in integrand                                  *
*   m (inp) - order of bessel function on centre a                     *
*   n (inp) - order of bessel function on centre b                     *
*   nt (inp) - number of integrals                                     *
*   maxord (inp) - maximum order of quadrature                         *
*   root (scr) - array of quadrature roots                             *
*   weight (scr) - array of quadrature weights                         *
*   work (scr) - work array for roots and weights evaluation           *
*   alpha (inp) - array of values a/2sqrt(c)                           *
*   beta (inp) - array of values b/2sqrt(c)                            *
*   gamma (inp) - array of values 1/sqrt(c)                            *
*   temp (scr) - scratch array                                         *
*   ind (scr) - index array                                            *
*   Q (out) - Q functions                                              *
*   tol (inp) - maximum relative error in bessel functions             *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer i,j,k,l,m,n,nt,ind(nt),order,maxord
      parameter (maxord = 12)
      double precision root(maxord),weight(maxord),work(maxord),
     &    temp(nt,12),alpha(nt),beta(nt),gamma(nt),Q(nt),r,zero,two,tol
      parameter (zero = 0.0d00, two = 2.0d00)
*
      if (nt .le. 0) return
      if (order.gt.maxord) call errquit(
     &      'ecp_t2_ghq: Order gt maxord diff is',(order-maxord),
     &       BASIS_ERR)
*
C     write (6,*) order
      call ecp_GAUSSQ (4,order,zero,zero,0,temp,work,root,weight)
Cbug
C      write (6,'(i5,1p2e20.10)') (i,root(i),weight(i),i=1,order)
Cend
      do i = 1,nt
        temp(i,1) = alpha(i)+beta(i)
        Q(i) = zero
      end do
*
      do l = 1,order
        j = 0
        do i = 1,nt
          r = root(l)+temp(i,1)
          if (r .ge. zero) then
            j = j+1
            ind(j) = i
            temp(j,2) = r**k
            temp(j,3) = alpha(i)*r*two
            temp(j,4) = beta(i)*r*two
          end if
        end do
        call ecp_bessel (j,m,temp(1,3),temp(1,5),temp(1,7),temp(1,6),
     &      tol)
        call ecp_bessel (j,n,temp(1,4),temp(1,6),temp(1,7),temp(1,3),
     &      tol)
        do i = 1,j
          Q(ind(i)) = Q(ind(i))+weight(l)*temp(i,2)*temp(i,5)*temp(i,6)
        end do
      end do
*
      do i = 1,nt
        Q(i) = Q(i)*gamma(i)**(k+1)
      end do
*
      return
      end
