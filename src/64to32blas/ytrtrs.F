      SUBROUTINE YTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
*
* $id$
*
      implicit none
#include "y64.fh"
      character          diag, trans, uplo
      integer            info, lda, ldb, n, nrhs
      double precision   a( lda, * ), b( ldb, * )
c
      INTRG4 n4,info4,lda4,lbd4,nrhs4
c
      n4=n
      nrhs4=nrhs
      lda4=lda
      ldb4=ldb
c
      call dtrtrs( uplo, trans, diag, n4, nrhs4, a, lda4, b, ldb4,
     $                   info4 )
      info=info4
      return
      end
