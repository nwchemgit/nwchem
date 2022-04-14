      subroutine ysyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
*
* $id$
*
      implicit none
#include "y64.fh"
      double precision alpha,beta
      integer k,lda,ldb,ldc,n
      character trans,uplo
      double precision a(lda,*),b(ldb,*),c(ldc,*)
c
      INTGR4 lda4,ldb4,ldc4
      INTGR4 n4,k4
c      
      lda4=lda
      ldb4=ldb
      ldc4=ldc
      n4=n
      k4=k
c      
      call dsyr2k(uplo,trans,n4,k4,alpha,a,lda4,b,ldb4,beta,c,ldc4)
      return
      end
