      subroutine ysysv( uplo, n, nrhs, a, lda, ipiv, b, ldb, work,
     $                  lwork, info )
*
* $Id$
*
      implicit none
#include "y64.fh"
      character          uplo
      integer            info, lda, ldb, lwork, n, nrhs
      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * ), work( * )

c
      INTGR4 n4,nrhs4,lda4,ldb4,lwork4,info4
c
      n4=n
      nrhs4=nrhs
      lda4=lda
      ldb4=ldb
      lwork4=lwork
      call dsysv( uplo, n4, nrhs4, a, lda4, ipiv, b, ldb4, work,
     $                  lwork4, info4 )
      info=info4
      return
      end
