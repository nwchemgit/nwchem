C $Id$
      Subroutine drdy_rotf (fmat,amat,scr,natom,ndim,maxdim)
C
C  Transform Hessian matrix from Gaussian standard orientation to
C    original orientation.
C
      Implicit none
      Integer natom,ndim,maxdim
      Double Precision fmat(maxdim,ndim),amat(3,3),scr(3)
C
      Integer i,j,jj,k,l
      Double Precision sum

      do i = 1,ndim
         do j = 1,natom
            jj = (j-1)*3
            do k = 1,3
               sum=0.0
               do l = 1,3
                  sum = sum + amat(k,l)*fmat(i,jj+l)
               enddo
               scr(k) = sum
            enddo
            do k = 1,3
               fmat(i,jj+k) = scr(k)
            enddo
         enddo
      enddo
      do i = 1,ndim
         do j = 1,natom
            jj = (j-1)*3
            do k = 1,3
               sum=0.0
               do l = 1,3
                  sum = sum + amat(k,l)*fmat(jj+l,i)
               enddo
               scr(k) = sum
            enddo
            do k = 1,3
               fmat(jj+k,i) = scr(k)
            enddo
         enddo
      enddo
      return
      end
