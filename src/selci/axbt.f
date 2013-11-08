C> \file axbt.f
C> Matrix multiply
C>
C> \ingroup selci
C> @{
C>
C> \brief Matrix multiply
C>
C> Performs the matrix-matrix multiplication
C> \f{eqnarray*}{
C>   C &=& A * B^T
C> \f}
C>
      subroutine selci_axbt(a,mrowa,b,mrowb,c,mrowc,ncol,nlink,nrow)
*
* $Id$
*
      implicit real*8 (a-h,o-z)
      integer mrowa !< [Input] The number of rows in A
      integer mrowb !< [Input] The number of rows in B
      integer mrowc !< [Input] The number of rows in C
      integer ncol  !< [Input] The length of a column in C
      integer nrow  !< [Input] The length of a row in C
      integer nlink !< [Input] The length of the contraction
      parameter (zero=0.0d0)
      double precision a(mrowa,nlink) !< [Input] Matrix A
      double precision b(mrowb,nlink) !< [Input] Matrix B
      double precision c(mrowc,nrow)  !< [Output] Matrix C
c     
c     matrix multiply c = a*bt
c     
*     mdc*if alliant
*     c(1:ncol,1:nrow) = matmul(a(1:ncol,1:nlink),
*     $     transpose(b(1:nrow,1:nlink)))
*     mdc*else
c     
c     Assumed that this is being used on small sparse matrices
c     
c     Optimize with loop unrolling
c     
      integer ind(3)
      real*8  bjk(3)
c
c     small cases
c
      if (nlink .eq. 1) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(j,1)
            enddo
         enddo
         return
      else if (nlink.eq.2) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(j,1) + a(i,2)*b(j,2)
            enddo
         enddo
         return
      else if (nlink .eq. 3) then
         do j = 1, nrow
            do i = 1, ncol
               c(i,j) = a(i,1)*b(j,1) + a(i,2)*b(j,2) + a(i,3)*b(j,3)
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
            test = b(j,k)
            if (test .ne. 0.0d0) then
               ndo = ndo + 1
               ind(ndo) = k
               bjk(ndo) = test
            endif
            if (ndo.eq.3) then
               k1 = ind(1)
               k2 = ind(2)
               k3 = ind(3)
               bjk1 = bjk(1)
               bjk2 = bjk(2)
               bjk3 = bjk(3)
               do i = 1,ncol
                  c(i,j)=c(i,j)+a(i,k1)*bjk1+a(i,k2)*bjk2+a(i,k3)*bjk3
               enddo
               ndo = 0
            endif
         enddo
         if (ndo.eq.2) then
            k1 = ind(1)
            k2 = ind(2)
            bjk1 = bjk(1)
            bjk2 = bjk(2)
            do i = 1,ncol
               c(i,j) = c(i,j) + a(i,k1)*bjk1 + a(i,k2)*bjk2
            enddo
         else if (ndo.eq.1) then
            k1 = ind(1)
            bjk1 = bjk(1)
            do i = 1,ncol
               c(i,j) = c(i,j) + a(i,k1)*bjk1
            enddo
         endif
c     
      enddo
*     *mdc*endif
*     c
*     call selci_mxma(a,1,mrowa,b,mrowb,1,c,1,mrowc,ncol,nlink,nrow)
*     
*     call dgemm( 'n', 't', ncol, nrow, nlink, 1.0d0, a, mrowa,
*     *     b, mrowb, 0.0d0, c, mrowc )
c     
      end
C>
C> @}
