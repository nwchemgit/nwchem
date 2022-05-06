      SUBROUTINE ygetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
*
* $id$
*
      implicit none
#include "y64.fh"
      character*1        trans
      integer            ipiv(*)
      integer            info, lda, ldb, n, nrhs, lwork
      double precision   a( lda, * ), b(ldb, * )
c
      INTGR4 n4, info4, lda4, ldb4, nrhs4
c
      n4 = n
      nrhs4 = nrhs
      lda4 = lda
      ldb4 = ldb
      lwork4 = lwork
c
      call dgetrs(trans, n4, nrhs4, a, lda4, ipiv, b, ldb4, info4)
      info=info4
      return
      end
