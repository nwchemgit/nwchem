C $Id: rotf.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_rotf (fmat,amat,scr,natom,ndim,ndima,maxdim)
C
C  Transform Hessian matrix from Gaussian standard orientation to
C    original orientation.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension fmat(maxdim,ndim),amat(ndima,ndima),scr(ndima)
C
      do i = 1,ndim
         do j = 1,natom
            jj = (j-1)*ndima
            do k = 1,ndima
               sum=0.0
               do l = 1,ndima
                  sum = sum + amat(k,l)*fmat(i,jj+l)
               enddo
               scr(k) = sum
            enddo
            do k = 1,ndima
               fmat(i,jj+k) = scr(k)
            enddo
         enddo
      enddo
      do i = 1,ndim
         do j = 1,natom
            jj = (j-1)*ndima
            do k = 1,ndima
               sum=0.0
               do l = 1,ndima
                  sum = sum + amat(k,l)*fmat(jj+l,i)
               enddo
               scr(k) = sum
            enddo
            do k = 1,ndima
               fmat(jj+k,i) = scr(k)
            enddo
         enddo
      enddo
      return
      end
            
