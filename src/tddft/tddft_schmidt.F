      subroutine tddft_schmidt(n,vector)
c
c $Id $
c
c Schmidt orthonormalize vectors
c Written by So Hirata, Feb 2002.
c (c) Battelle, PNNL, 2002.
c
      implicit none
      integer n
      double precision vector(n,n)
      integer i,j,k
      double precision a,b

      do i=1,n
        do j=i+1,n
          a=0.0d0
          do k=1,n
            a=a+vector(k,i)*vector(k,j)
          enddo
          b=0.0d0
          do k=1,n
            b=b+vector(k,i)*vector(k,i)
          enddo
          do k=1,n
            vector(k,j)=vector(k,j)-a/b*vector(k,i)
          enddo
        enddo
      enddo
      do i=1,n
        b=0.0d0
        do k=1,n
          b=b+vector(k,i)*vector(k,i)
        enddo
        do k=1,n
          vector(k,i)=vector(k,i)/dsqrt(b)
        enddo
      enddo
      return
      end
