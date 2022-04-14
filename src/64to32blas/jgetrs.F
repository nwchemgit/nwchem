      SUBROUTINE jgetrs( TRANS, N, NRHS, A, LDA, ipiv, B, LDB,
     $                   INFO )
*
* $id$
*
      implicit none
#include "y64.fh"
      character          trans
      integer ipiv(*)
      integer            info, lda, ldb, n, nrhs
      complex*16 a( lda, * ), b( ldb, * )
c
      INTGR4 n4,info4,lda4,ldb4,nrhs4
c
      n4=n
      nrhs4=nrhs
      lda4=lda
      ldb4=ldb
c
      call zgetrs( trans, n4, nrhs4, a, lda4, ipiv, b, ldb4,
     $                   info4 )
      info=info4
      return
      end
