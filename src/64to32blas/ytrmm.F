      subroutine ytrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
*
* $id$
*
      implicit none
#include "y64.fh"
      double precision alpha
      integer lda,ldb,m,n
      character diag,side,transa,uplo
      double precision a(lda,*),b(ldb,*)
c
      INTGR4 m4,n4
      INTGR4 lda4,ldb4
c      
      lda4=lda
      ldb4=ldb
      n4=n
      m4=m
c      
      
      call dtrmm(side,uplo,transa,diag,m4,n4,alpha,a,lda4,b,ldb4)
      return
      end
