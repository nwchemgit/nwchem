      subroutine ysytrs( uplo, n, nrhs, a, lda, ipiv, b, ldb, info )
*
* $id$
*
      implicit none
#include "y64.fh"
      character          uplo
      integer            info, lda, ldb, n, nrhs
      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )
c
      INTGR4 n4,info4,lda4,ldb4,nrhs4
c
      n4=n
      lda4=lda
      ldb4=ldb
      nrhs4=nrhs
c
      call dsytrs(uplo, n4, nrhs4, a, lda4, ipiv, b, ldb4, info4 )
      info=info4
      return
      end
