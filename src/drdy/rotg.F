C $Id$
      Subroutine drdy_rotg (g,amat,scr,natom,ndim)
C
C  Transform gradient vector from Gaussian standard orientation to
C    original orientation.
C
      Implicit none
      Integer natom,ndim
      Double Precision g(ndim),amat(3,3),scr(3)
C
      Integer i,ii,j,k
      Double Precision sum
C
      do i = 1,natom
         ii = 3*(i-1)
         do j = 1,3
            sum=0.0
            do k = 1,3
               sum = sum + amat(j,k)*g(ii+k)
            enddo
            scr(j) = sum
         enddo
         do j = 1,3
            g(ii+j) = scr(j)
         enddo
      enddo
      return
      end
