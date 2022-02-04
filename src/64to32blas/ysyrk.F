      subroutine ysyrk( uplo, trans, n, k, alpha, a, lda, beta, c, ldc )
*
* $Id$
*
      implicit none
#include "y64.fh"
      character            uplo, trans
      integer              lda, ldc, n, k
      double precision     alpha, beta, a( lda, * ), c( ldc, * )
c
      INTGR4 n4, k4, lda4, ldc4
c
      n4=n
      k4=k
      lda4=lda
      ldc4=ldc
      call dsyrk( uplo, trans, n4, k4, alpha, a, lda4, beta, c, ldc4 )
      return
      end

