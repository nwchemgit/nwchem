      subroutine ysymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
*
* $Id$
*
      implicit none
#include "y64.fh"
      double precision alpha,beta
      integer lda,ldb,ldc,m,n
      character side,uplo
*     ..
*     .. array arguments ..
      double precision a(lda,*),b(ldb,*),c(ldc,*)
c
      INTGR4 m4,n4,lda4,ldb4,ldc4
c
      m4=m
      n4=n
      lda4=lda
      ldb4=ldb
      ldc4=ldc
c
      call dsymm(side,uplo,m4,n4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
      return
      end
