C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_up_k (m_min,m_max,k,j,l,n,ldQ,alpha,beta,
     &    gamma,Q_lo,Q_mid,Q_hi)
*                                                                      *
*   Perform upward recursion in k for k even, using 6.38:              *
*   2c Q^{k+2}_{m,m+j} = (k-j) Q^k_{m,m+j} + a Q^{k+1}_{m+1,m+j}       *
*                    + b Q^{k+1}_{m,m-1+j}                             *
*   = gamma [(gamma (k-j)/2 Q^k_{m,m+j} + alpha Q^{k+1}_{m+1,m+j}      *
*                    + beta Q^{k+1}_{m,m-1+j} ]                        *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   m_min (inp) - minimum value of m                                   *
*   m_max (inp) - maximum value of m                                   *
*   k (inp) - value of k                                               *
*   j (inp) - value of j                                               *
*   l (inp) - switches alpha and beta for different forms of recursion *
*             formula.                                                 *
*   n (inp) - number of functions                                      *
*   alpha (inp) -  a/2*sqrt(c)                                         *
*   beta (inp) -  b/2*sqrt(c)                                          *
*   gamma (inp) - 1/sqrt(c)                                            *
*   Q_mid (inp) - Q^{k+1}_{m+1,m}                                      *
*   Q_lo (inp) - Q^{k}_{m,m}                                           *
*   Q_hi (out) - Q^{k+2}_{m,m}                                         *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,j,k,l,ldQ,m,m_min,m_max,n
      double precision two,fac,alpha(n),beta(n),gamma(n),
     &    Q_lo(ldQ,m_min:m_max),Q_mid(ldQ,m_min:m_max+1),
     &    Q_hi(ldQ,m_min:m_max)
*
      parameter (two = 2.0d00)
*
      fac = k-j
      fac = fac/two
C      write (6,'(/A,/)') 'upward_k_recursion'
C      write (6,*) 'k,j,m_min,m_max',k,j,m_min,m_max
      do m = m_min,m_max-j
C        write (6,*) 'm =',m
        do i = 1,n
          Q_hi(i,m) = gamma(i)*( gamma(i)*fac*Q_lo(i,m)
     &        + alpha(i)*Q_mid(i,m+l)+ beta(i)*Q_mid(i,m+1-l) )
C          write (6,'(1p4e20.10)') gamma(i),alpha(i),beta(i)
C          write (6,'(1p4e20.10)') Q_lo(i,m),Q_mid(i,m+l),Q_mid(i,m+1-l),
C     &        Q_hi(i,m)
        end do
      end do
*
      return
      end
