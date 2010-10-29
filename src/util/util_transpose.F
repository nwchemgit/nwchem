      subroutine util_transpose(a,lda,b,ldb,m,n)
*
* $Id$
*
      implicit none
c
      integer lda, ldb, m, n
      double precision a(lda,n), b(ldb,m)
c
c     Transpose a(m,n) to b(n,m)
c
c     The arrays CANNOT be the same.  Use util_inplace_transpose
c     instead.
c
      integer i, j
c$$$      integer b, ilo, jlo, ihi, jhi
c$$$      parameter (b=16)
c
      do j = 1, n
         do i = 1, m
            b(j,i) = a(i,j)
         enddo
      enddo
c$$$c
c$$$c     Untested blocking code
c$$$c
c$$$      do jlo = 1, n, b
c$$$         jhi = min(jlo+b-1,n)
c$$$         do ilo = 1, m, b
c$$$            ihi = min(ihi+b-1,m)
c$$$            do j = jlo, jhi
c$$$               do i = ilo, ihi
c$$$                  b(j,i) = a(i,j)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c
      end
      subroutine util_inplace_transpose(n, a)
      implicit none
c
      integer n
      double precision a(n,n)
c
      double precision aij, aji
c
      integer i, j
      do i = 1, n
         do j = 1, i-1
            aij = a(i,j)
            aji = a(j,i)
            a(i,j) = aji
            a(j,i) = aij
         enddo
      enddo
c
      end
      subroutine util_2_transpose(rows,cols,mat,ldmat,tmat,ldtmat)
      implicit none
#include "errquit.fh"
*
*
*  This routine forms the matrix transpose of mat.  The input 
*  matrix (mat) is a rows by cols matrix (stored in a ldmat 
*  by cols matrix in the calling routine).  The resultant matrix
*  (tmat) is a cols by rows matrix (stored in a ldtmat by rows 
*  matrix in the calling routine).  Matrices mat and tmat can
*  be the same storage but the routine must be called appropriately
*  (i.e., ldmat = ldtmat). 
*
*::passed:pointers
      integer rows   ! [input] number of rows in mat (columns in tmat)
      integer cols   ! [input] number of columns in mat (rows in tmat)
      integer ldmat  ! [input] leading dimension of mat
      integer ldtmat ! [input] leading dimension of tmat
      double precision   mat (ldmat ,cols) ! [input] matrix to be transposed
      double precision   tmat(ldtmat,rows) ! [output] transposed matrix
*::local:variables
* 
      integer r        ! loop index
      integer c        ! loop index
      integer square   ! square block size of non-square matrix
      double precision local_temp ! temporary variable for transpose
*
*:: do the square part of the matrix
      square = min(rows,cols)
      do r = 1,square
        do c = 1,r
          local_temp     = mat(c,r)
          tmat(c,r) = mat(r,c)
          tmat(r,c) = local_temp
        enddo
      enddo
      if      (rows.eq.cols) return ! square matrix return
*
      if      (rows.lt.cols) then   ! simple copy of non-square block
*
* e.g., 
*--------------------------------------
*   ***x    ***
*   ***y -> ***
*   ***z    ***
*           xyz
*--------------------------------------
*
        do r = (rows+1),cols
          do c = 1,rows
            tmat(r,c) = mat(c,r)
          enddo
        enddo
      else if (rows.gt.cols) then   ! simple copy of non-square block
*
* e.g., 
*--------------------------------------
*   ***    ***x
*   *** -> ***y
*   ***    ***z
*   xyz        
*--------------------------------------
*
        do r = 1,cols
          do c = (cols+1),rows
            tmat(r,c) = mat(c,r)
          enddo
        enddo
      else
        call errquit('util_transpose: square block did not return',911,
     &       UNKNOWN_ERR)
      endif
      end
