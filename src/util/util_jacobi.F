      subroutine util_jacobi(n, a, lda, e)
*
* $Id: util_jacobi.F,v 1.5 2003-10-17 22:58:35 carlfahl Exp $
*
      implicit none
      integer n, lda
      double precision a(lda,n)
      double precision e(n)
      call util_jacobi0(n, a, lda, e,.true.)
      return
      end
      subroutine util_jacobinosort(n, a, lda, e)
*
* $Id: util_jacobi.F,v 1.5 2003-10-17 22:58:35 carlfahl Exp $
*
      implicit none
      integer n, lda
      double precision a(lda,n)
      double precision e(n)
      call util_jacobi0(n, a, lda, e,.false.)
      return
      end
      subroutine util_jacobi0(n, a, lda, e,sorting)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer n, lda
      double precision a(lda,n)
      double precision e(n)
      logical sorting
c
c     Diagonalize A using Jacobi overwriting a with 
c     the eigenvectors and returning in e the eigenvalues
c
      integer l_v, k_v, i, j, ij, ierr
c
c     trivial 1x1 case
c
      if(n.eq.1) then
         e(1)=a(1,1)
         a(1,1)=1d0
         return
      endif
      if (.not. ma_push_get(mt_dbl, n*(n+1)/2 ,'v',l_v, k_v))
     $     call errquit('util_jacobi: v', n*(n+1)/2, MA_ERR)
c
c     Lower triangle of a into v
c
      ij = 0
      do i = 1, n
         do j = 1, i
            dbl_mb(k_v+ij) = a(j,i)
            ij = ij + 1
         enddo
      enddo
c
c     unit matrix into a
c
      call dfill(lda*n, 0d0, a, 1)
      call dfill(n, 1d0, a, lda+1)
      call stpr_sjacobi(lda, n, dbl_mb(k_v), e, a, ierr,sorting)
      if (ierr .ne. 0) call errquit('util_jacobi: ierr ', ierr,
     &       UNKNOWN_ERR)
c
      if (.not. ma_pop_stack(l_v)) call errquit('jacobi:ma',0, MA_ERR)
c
      end
