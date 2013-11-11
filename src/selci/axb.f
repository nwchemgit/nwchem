C> \ingroup selci
C> @{
C>
      subroutine selci_axb(a,mrowa,b,mrowb,c,mrowc,ncol,nlink,nrow)
*
* $Id$
*
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
      dimension a(mrowa,nlink),c(mrowc,nrow),b(mrowb,nrow)
c     
c     matrix multiply c = a*b
c     
c     Assumed that this is being used on small sparse matrices
c     
c     Optimize with loop unrolling
c     
      integer ind(3)
      real*8  bkj(3)
c
c     small cases
c
      if (nlink .eq. 1) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(1,j)
            enddo
         enddo
         return
      else if (nlink.eq.2) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j)
            enddo
         enddo
         return
      else if (nlink .eq. 3) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + a(i,3)*b(3,j)
            enddo
         enddo
         return
      endif
c
c     general case
c
      do j = 1,nrow
         do i = 1,ncol
            c(i,j) = zero
         enddo
c
         ndo = 0
         do k = 1,nlink
            test = b(k,j)
            if (test .ne. zero) then
               ndo = ndo + 1
               ind(ndo) = k
               bkj(ndo) = test
            endif
            if (ndo.eq.3) then
               k1 = ind(1)
               k2 = ind(2)
               k3 = ind(3)
               bkj1 = bkj(1)
               bkj2 = bkj(2)
               bkj3 = bkj(3)
               do i = 1,ncol
                  c(i,j)=c(i,j)+a(i,k1)*bkj1+a(i,k2)*bkj2+a(i,k3)*bkj3
               enddo
               ndo = 0
            endif
         enddo
         if (ndo.eq.2) then
            k1 = ind(1)
            k2 = ind(2)
            bkj1 = bkj(1)
            bkj2 = bkj(2)
            do i = 1,ncol
               c(i,j) = c(i,j) + a(i,k1)*bkj1 + a(i,k2)*bkj2
            enddo
         else if (ndo.eq.1) then
            k1 = ind(1)
            bkj1 = bkj(1)
            do i = 1,ncol
               c(i,j) = c(i,j) + a(i,k1)*bkj1
            enddo
         endif
c     
      enddo
c
      end
C>
C> @}
