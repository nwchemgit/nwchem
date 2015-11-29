C> \ingroup wfn1
C> @{
C>
C> \brief Calculate exp(K) where K is a matrix
C>
C> Calculate the exponent of a matrix by expanding the exponent
C> as a Taylor series. The algorithm uses only matrix-matrix
C> multiplications. 
C>
       subroutine wfn1_exp(n,K,m1,m2,expK)
       implicit none
c
       integer n !< [Input] The dimension of the matrices
       double precision K(n,n) !< [Input] The matrix K
       double precision m1(n,n) !< [Scratch] work space
       double precision m2(n,n) !< [Scratch] work space
       double precision expK(n,n) !< [Output] The matrix exponent of K
c
       double precision dnormK  !< The norm of K
       double precision dnormKn !< The norm of the n-th term in the
                                !< taylor series of K
       double precision dnorma  !< Alternative norm
c
       integer i, j, m, p, pi2
c
cDEBUG
c      write(*,*)'*** hvd: wfn1_exp: ',n
cDEBUG
c
c      Cap the rotations
c
       pi2 = 2.0d0*acos(-1.0d0)
       dnorma = 0.0d0
       do i = 1, n
         do j = 1, n
           dnorma = max(dnorma,abs(K(j,i)))
         enddo
       enddo
       if (dnorma.gt.pi2) then
         dnorma = pi2/dnorma
         do i = 1, n
           do j = 1, n
             K(j,i) = K(j,i)*dnorma
           enddo
         enddo
       endif
c
       dnormK = 0.0d0
       do i = 1, n
         do j = 1, n
           dnormK = dnormK + K(j,i)**2
         enddo
       enddo
       dnormK = sqrt(dnormK)
       dnormKn = sqrt(1.0d0*n)
       m = 0
       do while (dnormKn.gt.1.0d-20)
         m = m + 1
         dnormKn = dnormKn*dnormK/dble(m)
       enddo
c
       call dcopy(n*n,K,1,m1,1)
       call dscal(n*n,1.0d0/dble(m),m1,1)
       do i = 1, n
         m1(i,i) = m1(i,i) + 1.0d0
       enddo
       do p = m-1, 1, -1
         call dgemm('n','n',n,n,n,1.0d0/dble(p),K,n,m1,n,0.0d0,m2,n)
         call dcopy(n*n,m2,1,m1,1)
         do i = 1, n
           m1(i,i) = m1(i,i) + 1.0d0
         enddo
       enddo
       call dcopy(n*n,m1,1,expK,1)
cDEBUG
c      call dgemm('n','t',n,n,n,1.0d0,expK,n,expK,n,0.0d0,m2,n)
c      write(*,*)'*** hvd: wfn1_exp: unit mat'
c      call hess_hssout(m2,n,n,n)
cDEBUG
c
       end
C>
C> @}




