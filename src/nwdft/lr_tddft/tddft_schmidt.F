      subroutine tddft_schmidt(n,vector)
c
c $Id$
c
c Schmidt orthonormalize vectors
c Written by So Hirata, Feb 2002.
c (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
      integer n
      double precision vector(n,n)
      integer i,j,k
      double precision a

      do i=1,n
        a=0.0d0
        do k=1,n
          a=a+vector(k,i)*vector(k,i)
        enddo
        if (a.ne.0.d0) then
         do k=1,n
           vector(k,i)=vector(k,i)/dsqrt(a)
         enddo
        endif ! a ne 0.d0
        if (i.ne.n) then
          do j=i+1,n
            a=0.0d0
            do k=1,n
              a=a+vector(k,i)*vector(k,j)
            enddo
            do k=1,n
              vector(k,j)=vector(k,j)-a*vector(k,i)
            enddo
          enddo
        endif
      enddo
      return
      end
