C $Id: rotg.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_rotg (g,amat,scr,natom,ndim,ndima)
C
C  Transform gradient vector from Gaussian standard orientation to
C    original orientation.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension g(ndim),amat(ndima,ndima),scr(ndima)
C
      do i = 1,natom
         ii = ndima*(i-1)
         do j = 1,ndima
            sum=0.0
            do k = 1,ndima
               sum = sum + amat(j,k)*g(ii+k)
            enddo
            scr(j) = sum
         enddo
         do j = 1,ndima
            g(ii+j) = scr(j)
         enddo
      enddo
      return
      end
            
